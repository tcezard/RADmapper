'''
Created on 8 Apr 2010

@author: tcezard
'''
import MySQLdb
from utils import utils_param


####################################
#        database utilities        #
####################################

created_mysql_dbh={}
def getDbConnection(db='', user=None, pwd=None, host=None):
    """
    This is a general method to access a database.
    It get the database connection if it already exists or create it.
    """
    #mysql_dbh=created_mysql_dbh.get('%s%s%s'%(host,db,user))
    #if not mysql_dbh and user and pwd and host:  
    if user and pwd and host:  
        mysql_dbh = MySQLdb.connect(db = db,
                        host=host,
                        user=user,
                        passwd=pwd )
    return mysql_dbh

def getWtssDBConnection(db=''):
    """Get access to the database specified in the config file."""
    wtss_param=utils_param.get_wtss_parameter()
    mysql_dbh=getDbConnection(db, wtss_param.get_mysql_user(),
                              wtss_param.get_mysql_password(),
                              wtss_param.get_mysql_host())
    return mysql_dbh


def get_databases_name():
    from gene_annotation import annotation_db_def
    tables_def=annotation_db_def.get('tables')
    mysql_dbh=getWtssDBConnection()
    cursor=mysql_dbh.cursor()
    cursor.execute("show databases")
    res=cursor.fetchall()
    all_db=[]
    valid_db=[]
    for (db,) in res:
        if not db=="lost+found" and not db=="information_schema" and not db=="mysql" :
            all_db.append(db)
    for db in all_db:
        cursor.execute("use %s"%db)
        cursor.execute("show tables")
        res=cursor.fetchall()
        all_tables=[]
        for (table,) in res:
            all_tables.append(table)
        valid=True
        for table in tables_def:
            if not table in all_tables:
                valid=False
                break
        if valid:
            valid_db.append(db)
    return valid_db

