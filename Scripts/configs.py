# Import libraries
from configparser import ConfigParser


def config_db(filename, section='postgresql'):
    """
    config to connect to PostgreSQL server chemosimdb.
    """
    parser = ConfigParser()
    parser.read(filename)

    res = {}
    if parser.has_section(section):
        params = parser.items(section)
        res.update({param[0] : param[1] for param in params})
    else:
        raise Exception('Section {0} not found in the {1} file'.format(section, filename))
    return res


# if __name__ =='__main__':
#     filename = ''
#     result = config_db()
#     print(result)