from .main import final_iteration, filter_results
from .utils import prot_gen
import pandas as pd
from collections import defaultdict, Counter
import argparse
import logging
import ete3
from ete3 import NCBITaxa
ncbi = NCBITaxa()

logger = logging.getLogger(__name__)

def run():
    parser = argparse.ArgumentParser(
        description='Combine DirectMS1 search results',
        epilog='''

    Example usage
    -------------
    $ ms1combine.py file1_PFMs_ML.tsv ... filen_PFMs_ML.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', nargs='+', help='input tsv PFMs_ML files for union')
    parser.add_argument('-d', '-db', help='path to protein fasta file', required=True)
    parser.add_argument('-out', help='prefix output file names', default='group_specific_statistics_by_')
    parser.add_argument('-prots_full', help='path to any of *_proteins_full.tsv file. By default this file will be searched in the folder with PFMs_ML files', default='')
    parser.add_argument('-fdr', help='protein fdr filter in %%', default=1.0, type=float)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-nproc', help='number of processes', default=1, type=int)
    parser.add_argument('-groups', help="dbname: To use taxonomy in protein name. OX: Use OX= from fasta file. Or can be 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain'", default='dbname')
    parser.add_argument('-pp', help='protein priority table for keeping protein groups when merge results by scoring', default='')
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)

    group_to_use = args['groups']
    allowed_groups = [
        'dbname', 'OX', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain'
    ]
    if group_to_use not in allowed_groups:
        logging.critical('group is not correct! Must be: %s', ','.join(allowed_groups))
        return -1

    dbname_map = dict()
    ox_map = dict()
    for dbinfo, dbseq in prot_gen(args):
        dbname = dbinfo.split(' ')[0]

        if group_to_use != 'dbname':
            try:
                ox = dbinfo.split('OX=')[-1].split(' ')[0]
            except:
                ox = 'Unknown'
        else:
            try:
                ox = dbinfo.split(' ')[0].split('|')[-1].split('_')[-1]
            except:
                ox = 'Unknown'
        dbname_map[dbname] = ox

    cnt = Counter(dbname_map.values())

    if group_to_use not in ['dbname', 'OX']:
        for ox in cnt.keys():

            line = ncbi.get_lineage(ox)
            ranks = ncbi.get_rank(line)
            if group_to_use not in ranks.values():
                logger.warning('%s does not have %s', str(ox), group_to_use)
                group_custom = 'OX:' + ox
                # print('{} does not have {}'.format(i, group_to_use))
                # continue

            else:
                ranks_rev = {k[1]:k[0] for k in ranks.items()}
                # print(ranks_rev)
                group_custom = ranks_rev[group_to_use]

            ox_map[ox] = group_custom


        for dbname in list(dbname_map.keys()):
            dbname_map[dbname] = ox_map[dbname_map[dbname]]

        cnt = Counter(dbname_map.values())

    print(cnt.most_common())

    # return -1

    d_tmp = dict()

    df1 = None
    for idx, filen in enumerate(args['file']):
        logging.info('Reading file %s' % (filen, ))
        df3 = pd.read_csv(filen, sep='\t', usecols=['ids', 'qpreds', 'preds', 'decoy', 'seqs', 'proteins', 'peptide', 'iorig'])
        df3['ids'] = df3['ids'].apply(lambda x: '%d:%s' % (idx, str(x)))
        df3['fidx'] = idx

        df3 = df3[df3['qpreds'] <= 10]


        qval_ok = 0
        for qval_cur in range(10):
            if qval_cur != 10:
                df1ut = df3[df3['qpreds'] == qval_cur]
                decoy_ratio = df1ut['decoy'].sum() / len(df1ut)
                d_tmp[(idx, qval_cur)] = decoy_ratio
                # print(filen, qval_cur, decoy_ratio)

        if df1 is None:
            df1 = df3
            if args['prots_full']:
                df2 = pd.read_csv(args['prots_full'], sep='\t')
            else:
                try:
                    df2 = pd.read_csv(filen.replace('_PFMs_ML.tsv', '_proteins_full.tsv'), sep='\t')
                except:
                    logging.critical('Proteins_full file is missing!')
                    break

        else:
            df1 = pd.concat([df1, df3], ignore_index=True)

    d_tmp = [z[0] for z in sorted(d_tmp.items(), key=lambda x: x[1])]
    qdict = {}
    for idx, val in enumerate(d_tmp):
        qdict[val] = int(idx / len(args['file']))
    df1['qpreds'] = df1.apply(lambda x: qdict[(x['fidx'], x['qpreds'])], axis=1)


    pept_prot = defaultdict(set)
    group_to_pep = defaultdict(set)
    for seq, prots in df1[['seqs', 'proteins']].values:
        for dbname in prots.split(';'):
            pept_prot[seq].add(dbname)
            group_to_pep[dbname_map[dbname]].add(seq)

    protsN = dict()
    for dbname, theorpept in df2[['dbname', 'theoretical peptides']].values:
        protsN[dbname] = theorpept


    prefix = args['prefix']
    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]
    fdr = float(args['fdr']) / 100

    # all_proteins = []

    base_out_name = args['out'] + group_to_use + '.tsv'

    out_dict = dict()

    for group_name in cnt:

        logging.info(group_name)

        df1_tmp = df1[df1['peptide'].apply(lambda x: x in group_to_pep[group_name])]

        protsN_tmp = dict()
        for k, v in protsN.items():
            dbname_map[k] == group_name
            protsN_tmp[k] = v

        resdict = dict()

        resdict['qpreds'] = df1_tmp['qpreds'].values
        resdict['preds'] = df1_tmp['preds'].values
        resdict['seqs'] = df1_tmp['peptide'].values
        resdict['ids'] = df1_tmp['ids'].values
        resdict['iorig'] = df1_tmp['iorig'].values

        # mass_diff = resdict['qpreds']
        # rt_diff = resdict['qpreds']

        e_ind = resdict['qpreds'] <= 9
        resdict = filter_results(resdict, e_ind)

        mass_diff = resdict['qpreds']
        rt_diff = resdict['qpreds']

        if args['pp']:
            df4 = pd.read_table(args['pp'])
            prots_spc_basic2 = df4.set_index('dbname')['score'].to_dict()
        else:
            prots_spc_basic2 = False


    
        top_proteins = final_iteration(resdict, mass_diff, rt_diff, pept_prot, protsN_tmp, base_out_name, prefix, isdecoy, isdecoy_key, escore, fdr, args['nproc'], prots_spc_basic2=prots_spc_basic2, output_all=False)
        # all_proteins.extend(top_proteins)
        out_dict[group_name] = len(top_proteins)
        # print(top_proteins)
        print('\n')
        # break
    
    with open(base_out_name, 'w') as output:
        output.write('taxid\tproteins\n')
        for k, v in out_dict.items():
            output.write('\t'.join((str(k), str(v))) + '\n')

    

if __name__ == '__main__':
    run()
