import argparse
import os
import numpy as np
from scipy import stats
import xlrd, xlwt
from xlutils.copy import copy

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('inhibitor', type=str)
    parser.add_argument('-sites', type=str, nargs='*', required=True)
    return parser.parse_args()

def calculate_Z_score(sum_workbook, inhibitor, sites, row):
    sum_sheet = sum_workbook.get_sheet(inhibitor)
    for site in sites:
        workbook_name = inhibitor + '/' + site + '.xls'
        workbook_read = xlrd.open_workbook(workbook_name)
        for aa in ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', \
                'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']:
            mutation = site + aa
            sheet_read = workbook_read.sheet_by_name(mutation)
            ddG_fold = sheet_read.cell_value(1, 1)
            ddG_bind_total = np.array(sheet_read.col_values(2)[2:]).astype(np.float64)
            ddG_substrate = np.array(sheet_read.col_values(5)[2:]).astype(np.float64)
            ddG_bind_total_zscores = stats.zscore(ddG_bind_total)
            # first_max_ddG_bind_total, second_max_ddG_bind_total = sorted(ddG_bind_total[:-1], reverse=True)[:2]
            # if ddG_bind_total[-1] > max_ddG_bind_total and max_ddG_bind_total < 5 and ddG_fold < 15:
            if True: #ddG_bind_total_zscores[-1] > 0.75:
                sum_sheet.write(row, 0, mutation)
                sum_sheet.write(row, 1, ddG_fold)
                sum_sheet.write(row, 2, ddG_bind_total[-1])
                sum_sheet.write(row, 3, ddG_substrate[-1])
                sum_sheet.write(row, 4, ddG_bind_total_zscores[-1])
                for col in range(len(ddG_bind_total) - 1):
                    sum_sheet.write(row, col * 2 + 5, ddG_bind_total[col])
                    sum_sheet.write(row, col * 2 + 6, ddG_substrate[col])
                    # sum_sheet.write(row, col * 2 + 5, ddG_bind_total_zscores[col])
                row += 1

if __name__ == '__main__':
    args = parse_arguments()
    workbook_name = args.inhibitor + '/' + args.inhibitor + '.xls'
    if not os.path.isfile(workbook_name):
        sum_workbook = xlwt.Workbook(encoding="ascii")
        sum_sheet = sum_workbook.add_sheet(args.inhibitor)
        sum_sheet.write(0, 0, 'mutation')
        sum_sheet.write(0, 1, 'ddG_fold')
        sum_sheet.write(0, 2, args.inhibitor + '_ddG_bind')
        sum_sheet.write(0, 3, args.inhibitor + '_ddG_sub')
        sum_sheet.write(0, 4, 'z_score')
        col = 5
        for sub in filter(lambda x: x != 11, range(4, 16)):
            if sub in [7, 12, 13]:
                sum_sheet.write(0, col, str(sub) + '-' + str(sub + 1) + '_ddG_bind')
                sum_sheet.write(0, col + 1, str(sub) + '-' + str(sub + 1) + '_ddG_sub')
                col += 2
            else:
                sum_sheet.write(0, col, str(sub) + '-' + str(sub + 1) + '_ddG_bind')
                sum_sheet.write(0, col + 1, str(sub) + '-' + str(sub + 1) + '_ddG_sub')
                sum_sheet.write(0, col + 2, str(sub) + '-' + str(sub + 1) + '_ddG_bind')
                sum_sheet.write(0, col + 3, str(sub) + '-' + str(sub + 1) + '_ddG_sub')
                col += 4
        row = 1
    else:
        sum_workbook_read = xlrd.open_workbook(workbook_name)
        sum_sheet_read = sum_workbook_read.sheet_by_index(0)
        row = sum_sheet_read.nrows
        sum_workbook = copy(sum_workbook_read)
    calculate_Z_score(sum_workbook, args.inhibitor, args.sites, row)
    sum_workbook.save(workbook_name)
