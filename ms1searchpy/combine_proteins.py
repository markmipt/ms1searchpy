from pyteomics import auxiliary as aux
import pandas as pd
import argparse
import logging

logger = logging.getLogger(__name__)

def run():
    parser = argparse.ArgumentParser(
        description='Combine DirectMS1 search results',
        epilog='''

    Example usage
    -------------
    $ ms1combine_proteins file1_proteins_full.tsv ... filen_proteins_full.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', nargs='+', help='input tsv proteins_full files')
    parser.add_argument('-out', help='prefix for joint file name', default='combined')
    parser.add_argument('-fdr', help='protein fdr filter in %%', default=1.0, type=float)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO)

    tmp_list = []
    for idx, filen in enumerate(args['file']):
        tmp_list.append(pd.read_csv(filen, sep='\t'))
    
    df1 = pd.concat(tmp_list).groupby(['dbname']).sum().reset_index()
    out_name = '%s.features_proteins.tsv' % (args['out'], )




    prots_spc = df1.set_index('dbname').score.to_dict()
    fdr = args['fdr'] / 100

    prefix = args['prefix']
    isdecoy = lambda x: x[0].startswith(prefix)
    isdecoy_key = lambda x: x.startswith(prefix)
    escore = lambda x: -x[1]



    checked = set()
    for k, v in list(prots_spc.items()):
        if k not in checked:
            if isdecoy_key(k):
                if prots_spc.get(k.replace(prefix, ''), -1e6) > v:
                    del prots_spc[k]
                    checked.add(k.replace(prefix, ''))
            else:
                if prots_spc.get(prefix + k, -1e6) > v:
                    del prots_spc[k]
                    checked.add(prefix + k)

    filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True, correction=1)
    if len(filtered_prots) < 1:
        filtered_prots = aux.filter(prots_spc.items(), fdr=fdr, key=escore, is_decoy=isdecoy, remove_decoy=True, formula=1, full_output=True, correction=0)
    identified_proteins = 0

    for x in filtered_prots:
        identified_proteins += 1

    logger.info('TOP 5 identified proteins:')
    logger.info('dbname\tscore')
    for x in filtered_prots[:5]:
        logger.info('\t'.join((str(x[0]), str(x[1]))))
    logger.info('Final joint search: identified proteins = %d', len(filtered_prots))

    df1 = pd.DataFrame.from_records(filtered_prots, columns=['dbname', 'score'])
    df1.to_csv(out_name, index=False, sep='\t')

if __name__ == '__main__':
    run()
