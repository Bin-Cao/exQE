"""
exQE Package

This package provides tools for processing crystal structure in QE.
"""
import datetime
# Version of the exQE package
__version__ = '0.0.1'
__email__ = 'binjacobcao@gmail.com'

now = datetime.datetime.now()
formatted_date_time = now.strftime('%Y-%m-%d %H:%M:%S')
print('exQE, Bin CAO, HKUST.GZ, http://www.caobin.asia/' )
print('Executed on :',formatted_date_time, ' | Have a great day.')  
print('='*80)