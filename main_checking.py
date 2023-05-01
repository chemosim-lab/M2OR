import argparse
import os
import pandas
import datetime

from formatting import PreFormatter, PostFormatter
from checking import Checker, PostChecker, OptionalChecker
    


def main_check(csv_path, sep = ';', run_optional_checker = True, auxillary_dir = 'Data', log_dir = 'logs'):
    """
    main script to run checks and format the data.
    
    Parameters:
    -----------
    csv_path : str
        path to the csv to check.

    sep : str
        csv separator.

    run_optional_checker : bool
        whether to run optional checker. See OptionalChecker in checking.py for more details.

    Return:
    -------
    df : pandas.DataFrame
        formated dataframe. If the checks are not passed an error is raised.
    """
    df = pandas.read_csv(csv_path, sep = sep, index_col = 0)
    
    formatter = PreFormatter(auxillary_dir = auxillary_dir, log_dir = log_dir)
    df, exclude_df = formatter(df)
    
    # exclude_df.to_csv('RawData_test/exclude.csv', sep=';')
    checker = Checker(auxillary_dir = auxillary_dir, log_dir = log_dir)
    checker(df)
    
    post_formatter = PostFormatter(auxillary_dir = auxillary_dir, log_dir = log_dir)
    df = post_formatter(df)
    
    post_checker = PostChecker(auxillary_dir = auxillary_dir, log_dir = log_dir)
    post_checker(df)

    if run_optional_checker:
        optional_checker = OptionalChecker(auxillary_dir = auxillary_dir, log_dir = log_dir)
        optional_checker(df)

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_path', type=str, required=True,
                        help='path to csv for batch upload.')
    parser.add_argument('--sep', type=str, default=';',
                        help='separator for the csv. Semicolon is used by default')
    parser.add_argument('--additional_check', type=str, default='y',
                        help='whether to run additional check or not. y/n. It is run by deafult.')
    args = parser.parse_args()

    print('csv path: {}'.format(args.csv_path))
    
    csv_path = args.csv_path
    sep = args.sep
    if args.additional_check.lower() == 'y' or args.additional_check.lower() == 'yes':
        run_optional_checker = True
    else:
        run_optional_checker = False

    df = main_check(csv_path, sep, run_optional_checker)