import os
import sys
import argparse
import re
from datetime import datetime
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(description="Fetch the latest version of UKB data field.")
    parser.add_argument("--fids", nargs="+", help="UKB data field ids.")
    parser.add_argument("--baskets-dir", default="/tsd/p33/data/durable/s3-api/ukbdata/phenotypes/Baskets", help="Directory with UKB baskets.")
    parser.add_argument("--verbose", action='store_true', help="Print extended information.")
    parser.add_argument("--out", help="Output file. If not specified noting is extracted and saved.")
    return parser.parse_args(args)


class FieldID(object):
    """
    An object to handle UKB field ids.
    Field ID must have format:
        fid[-inst[.ind]]
    where fid/inst/ind are arbitrary integers. Parts in the square brakets are
    optional.
    """
    def __init__(self, id_str):
        fid_match = re.match(r'(?P<fid>\d+)', id_str)
        inst_match = re.match(r'\d+-(?P<inst>\d+)', id_str)
        ind_match = re.match(r'\d+-\d+\.(?P<ind>\d+)', id_str)
        if fid_match:
            self.fid = fid_match['fid']
            self.inst = inst_match['inst'] if inst_match else None
            self.ind = ind_match['ind'] if ind_match else None
        else:
            raise ValueError(f'Incorrect Field ID string: {id_str}')

    def __repr__(self):
        return f'{self.fid}-{self.inst if self.inst else ""}.{self.ind if self.ind else ""}'

    def __eq__(self, other):
        return ( (self.fid == other.fid) and
                 (self.inst is None or other.inst is None or self.inst == other.inst) and
                 (self.ind is None or other.ind is None or self.ind == other.ind) )

    
class BasketReleaseFile(object):
    def __init__(self, csv_path, basket_id, release_date, fids):
        self.csv_path = csv_path
        self.basket_id = basket_id
        self.release_date = release_date
        self.fids = fids
    
    @classmethod
    def from_prefix(cls, prefix):
        """
        Args:
            prefix: file path prefix for .csv and .log files produced by
                    ./ukb_conv *.enc_ukb csv
        """
        csv_path = f'{prefix}.csv'
        log_path = f'{prefix}.log'
        if os.path.isfile(csv_path) and os.path.isfile(log_path):
            basket_pattern = re.compile(r'Basket (?P<basket_id>\d+)')
            dt_pattern = re.compile(r'\d{4}-\d{2}-\d{2}.{1}\d{2}:\d{2}:\d{2}')
            basket_id = None
            release_date = None
            with open(log_path) as log_file:
                for l in log_file:
                    basket_match = basket_pattern.search(l)
                    dt_match = dt_pattern.search(l)
                    if basket_match and basket_id is None:
                        basket_id = basket_match['basket_id']
                    if dt_match and release_date is None:
                        release_date = datetime.fromisoformat(dt_match.group(0))
                    if basket_id and release_date:
                        break
            with open(csv_path) as csv_file:
                fids = [FieldID(fid_str) for fid_str in csv_file.readline().replace('"','').split(',')[1:]]
        else:
            raise ValueError(f'{csv_path} and/or {log_path} files do not exist.')
        return cls(csv_path, basket_id, release_date, fids)
    
    def __repr__(self):
        return f'Path: {self.csv_path}; Basket: {self.basket_id}; Release: {self.release_date}; #Fids = {len(self.fids)}; Fids: {self.fids}'
    
    
def get_basket_release_files(baskets_dir):
    print(f"Scanning baskets directory: {baskets_dir}")
    ukb_csv_pattern = re.compile(r'ukb\d+\.csv')
    # take only those release files which have non-empty csv files
    ukb_prefixes = [os.path.splitext(de)[0] for de in os.scandir(baskets_dir)
                    if ukb_csv_pattern.match(de.name) and de.stat().st_size > 0]
    print(f'    {len(ukb_prefixes)} basket release files')
    basket_release_files = [BasketReleaseFile.from_prefix(prefix) for prefix in ukb_prefixes]
    return basket_release_files


def get_basket_release_files_for_fid(fid, basket_release_files):
    """
    Args:
        fid:                  string or FieldID instance
        basket_release_files: list of BasketReleaseFile instances
    Return:
        list of BasketReleaseFile instances containing fid
    """
    if not isinstance(fid, FieldID):
        fid = FieldID(fid)
    basket_release_files_for_fid = []
    for brf in basket_release_files:
        relevant_fids = [i for i in brf.fids if i == fid]
        if relevant_fids:
            brf_for_fid = BasketReleaseFile(brf.csv_path, brf.basket_id, brf.release_date, relevant_fids)
            basket_release_files_for_fid.append(brf_for_fid)
    print(f'{len(basket_release_files_for_fid)} basket release files containing fid {fid}')
    return basket_release_files_for_fid



if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    basket_release_files = get_basket_release_files(args.baskets_dir)

    out_dfs = []
    fetched_fids = []
    for fid in args.fids:
        print(f'Processing {fid}')
        basket_release_files_for_fid = get_basket_release_files_for_fid(fid, basket_release_files)
        if args.verbose:
            print('\n'.join([f'    {i+1}: {brf}' for i,brf in enumerate(basket_release_files_for_fid)]))

        latest_brf_for_fid = max(basket_release_files_for_fid, key=lambda brf: brf.release_date)
        print(f'    Latest release file:\n        {latest_brf_for_fid}')

        if args.out:
            # extract data
            current_fids = []
            for fid_str in map(str,latest_brf_for_fid.fids):
                if fid_str in fetched_fids:
                    print(f'        {fid_str} has been fetched already, skipping now.')
                else:
                    current_fids.append(fid_str)
            if current_fids:
                fetched_fids += current_fids
                usecols = ['eid'] + current_fids
                df = pd.read_csv(latest_brf_for_fid.csv_path, usecols=usecols, index_col='eid', dtype=object, na_filter=False)
                out_dfs.append(df)

    if args.out:
        # merge and save
        df = out_dfs[0].join(out_dfs[1:], how='outer')
        df.to_csv(args.out)
        print(f'Merged results saved to: {args.out}')



