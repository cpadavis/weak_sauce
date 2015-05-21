"""

"""
from subprocess import call
from os import path, makedirs
import pandas as pd


def check_make(folder):
    if not path.exists(folder):
        makedirs(folder)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--job',
                        dest='job',
                        #required=True,
                        default='paths',
                        help='[paths, images]')
    parser.add_argument('-o',
                        dest='out',
                        default='/nfs/slac/g/ki/ki18/des/cpd/flats')
    parser.add_argument('--user',
                        dest='user')
    parser.add_argument('--password',
                        dest='password')
    options = parser.parse_args()
    user = options.user
    password = options.password
    out_directory = options.out
    csv_path = '{0}/httppath_exposure.csv'.format(out_directory)
    if options.job == 'paths':
        print('Downloading CSV of exposures.')
        query = "select path, obstype, band, exposurename, id, expnum, nite, mjd_obs from httppath_exposure where obstype in ('dome flat', 'sky flat')"

        try:
            import easyaccess
            connection = easyaccess.connect()
            cursor = connection.cursor()
            #csv = connection.query_to_pandas(query)
            connection.query_and_save(query, csv_path)
            connection.close()

        except:
            # huh. Well let's try commandline
            command = ['easyaccess', '-c',
                       "{0}; > {1}".format(query, csv_path)]
            call(command)

    elif options.job == 'images':
        print('Using CSV at {0}'.format(csv_path))
        csv = pd.read_csv(csv_path)
        # download based on band and obstype and cutoff date

        for obstype in csv['OBSTYPE'].unique():
            for band in csv['BAND'].unique():
                conds = (csv['OBSTYPE'] == obstype) * (csv['BAND'] == band)
                for index, row in csv[conds].iterrows():
                    httppath = row['PATH']
                    exposurename = row['EXPOSURENAME']

                    out_dir = '{0}/{1}/{2}'.format(out_directory, obstype, band)
                    check_make(out_dir)
                    command = ['bsub', '-W', str(5),
                               'wget',
                               '--user', user,
                               '--password', password,
                               '--no-check-certificate',
                               '-O', '{0}/{1}.fz'.format(out_dir, exposurename),
                               httppath]
                    call(command)
    else:
        print('Unrecognized job {0}! Exiting'.format(options.job))
