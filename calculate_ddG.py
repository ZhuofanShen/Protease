import argparse
import os
import numpy as np
import xlrd, xlwt
from xlutils.copy import copy


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('site_aa', type=str, nargs=2)
    parser.add_argument('-inh', '--inhibitor', type=str)
    return parser.parse_args()

def read_variant_delta_G(site, native_aa):
    # Read wild type apo protease delta G
    delta_G_data_file = 'apo_6YB7/' + site + '/' + site + '_'+ native_aa + '.dat'
    if os.path.isfile(delta_G_data_file):
        with open(delta_G_data_file) as pf:
            lines = pf.readlines()
        if len(lines) == 1:
            wt_apo_delta_G_data = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
        else:
            wt_apo_delta_G_data = (np.array(lines[0][:-1].split(',')).astype(np.float64), \
                    np.array(lines[1][:-1].split(',')).astype(np.float64))
    else:
        wt_apo_delta_G_data = None
    # Read wild type protease-nsp complex delta G and wild type apo protease delta G
    wt_delta_G_dict = dict()
    wt_apo_delta_G_dict = dict()
    for nsp in ['nsp4-nsp5_7T70', 'nsp5-nsp6_7T8M', 'nsp6-nsp7_7MB6', 'nsp7-nsp8_7T8R', \
            'nsp8-nsp9_7T9Y', 'nsp9-nsp10_7TA4', 'nsp10-nsp11_7TA7', 'nsp12-nsp13_7TB2', \
            'nsp13-nsp14_7TBT', 'nsp14-nsp15_7TA4', 'nsp15-nsp16_7TC4']:
        delta_G_data_file = nsp + '/' + site + '/' + site + '_'+ native_aa + '.dat'
        if os.path.isfile(delta_G_data_file):
            with open(delta_G_data_file) as pf:
                lines = pf.readlines()
            if len(lines) == 1:
                wt_delta_G_dict[nsp] = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
            else:
                wt_delta_G_dict[nsp] = (np.array(lines[0][:-1].split(',')).astype(np.float64), \
                        np.array(lines[1][:-1].split(',')).astype(np.float64))
        apo_state = 'apo_' + nsp.split('_')[1]
        apo_delta_G_data_file = apo_state + '/' + site + '/' + site + '_'+ native_aa + '.dat'
        if os.path.isfile(apo_delta_G_data_file):
            with open(apo_delta_G_data_file) as pf:
                lines = pf.readlines()
            if len(lines) == 1:
                wt_apo_delta_G_dict[apo_state] = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
            else:
                wt_apo_delta_G_dict[apo_state] = (np.array(lines[0][:-1].split(',')).astype(np.float64), \
                        np.array(lines[1][:-1].split(',')).astype(np.float64))
    return wt_apo_delta_G_data, wt_apo_delta_G_dict, wt_delta_G_dict

def calculate_variant_ddG(site, wt_apo_delta_G_data, wt_apo_delta_G_dict, wt_delta_G_dict):
    workbook = xlwt.Workbook(encoding="ascii")
    # for each sheet
    for aa in ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', \
                'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']:
        mutation = site + aa
        try:
            sheet = workbook.get_sheet(mutation)
        except:
            sheet = workbook.add_sheet(mutation, cell_overwrite_ok=True)
        sheet.write(0, 1, 'ddG_total')
        sheet.write(0, 2, 'ddG_bind,total')
        sheet.write(0, 3, 'ddG_shell')
        sheet.write(0, 4, 'ddG_bind,shell')
        sheet.write(0, 5, 'ddG_sub')
        sheet.write(0, 6, 'ddG_res')
        sheet.write(0, 7, 'ddG_bind,res')
        sheet.write(0, 8, 'ddG_interact')
        sheet.write(0, 9, 'ddG_cst')
        sheet.write(1, 0, 'apo')
        # Calculate apo protease ddG for each amino acid type
        apo_delta_G_data_file = 'apo_6YB7/' + site + '/' + site + '_'+ aa + '.dat'
        if os.path.isfile(apo_delta_G_data_file):
            with open(apo_delta_G_data_file) as pf:
                lines = pf.readlines()
            if len(lines) == 1:
                apo_delta_G_data = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
                apo_ddG_data = np.subtract(apo_delta_G_data, wt_apo_delta_G_data)
                sheet.write(1, 1, apo_ddG_data[0])
                sheet.write(1, 2, 'NA')
                sheet.write(1, 3, apo_ddG_data[1])
                sheet.write(1, 4, 'NA')
                sheet.write(1, 5, apo_ddG_data[2])
                sheet.write(1, 6, apo_ddG_data[3])
                sheet.write(1, 7, 'NA')
                sheet.write(1, 8, apo_ddG_data[4])
                sheet.write(1, 9, apo_ddG_data[5])
                row = 2
            else:
                apo_delta_G_data_1 = np.array(lines[0][:-1].split(',')).astype(np.float64)
                apo_ddG_data_1 = np.subtract(apo_delta_G_data_1, wt_apo_delta_G_data[0])
                sheet.write(1, 1, apo_ddG_data_1[0])
                sheet.write(1, 2, 'NA')
                sheet.write(1, 3, apo_ddG_data_1[1])
                sheet.write(1, 4, 'NA')
                sheet.write(1, 5, apo_ddG_data_1[2])
                sheet.write(1, 6, apo_ddG_data_1[3])
                sheet.write(1, 7, 'NA')
                sheet.write(1, 8, apo_ddG_data_1[4])
                sheet.write(1, 9, apo_ddG_data_1[5])
                apo_delta_G_data_2 = np.array(lines[1][:-1].split(',')).astype(np.float64)
                apo_ddG_data_2 = np.subtract(apo_delta_G_data_2, wt_apo_delta_G_data[1])
                sheet.write(2, 1, apo_ddG_data_2[0])
                sheet.write(2, 2, 'NA')
                sheet.write(2, 3, apo_ddG_data_2[1])
                sheet.write(2, 4, 'NA')
                sheet.write(2, 5, apo_ddG_data_2[2])
                sheet.write(2, 6, apo_ddG_data_2[3])
                sheet.write(2, 7, 'NA')
                sheet.write(2, 8, apo_ddG_data_2[4])
                sheet.write(2, 9, apo_ddG_data_2[5])
                row = 3
        else:
            row = 1
        # Calculate holo protease ddG for each amino acid type
        for nsp in ['nsp4-nsp5_7T70', 'nsp5-nsp6_7T8M', 'nsp6-nsp7_7MB6', 'nsp7-nsp8_7T8R', \
                'nsp8-nsp9_7T9Y', 'nsp9-nsp10_7TA4', 'nsp10-nsp11_7TA7', 'nsp12-nsp13_7TB2', \
                'nsp13-nsp14_7TBT', 'nsp14-nsp15_7TA4', 'nsp15-nsp16_7TC4']:
            delta_G_data_file = nsp + '/' + site + '/' + site + '_'+ aa + '.dat'
            apo_state = 'apo_' + nsp.split('_')[1]
            apo_delta_G_data_file = apo_state + '/' + site + '/' + site + '_'+ aa + '.dat'
            calculate_delta_G_bind = False
            if os.path.isfile(apo_delta_G_data_file):
                calculate_delta_G_bind = True
                with open(apo_delta_G_data_file) as pf:
                    apo_lines = pf.readlines()
            if os.path.isfile(delta_G_data_file):
                with open(delta_G_data_file) as pf:
                    lines = pf.readlines()
                if len(lines) == 1:
                    sheet.write(row, 0, nsp)
                    delta_G_data = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
                    ddG_data = np.subtract(delta_G_data, wt_delta_G_dict.get(nsp))
                    sheet.write(row, 1, ddG_data[0])
                    sheet.write(row, 3, ddG_data[1])
                    sheet.write(row, 5, ddG_data[2])
                    sheet.write(row, 6, ddG_data[3])
                    sheet.write(row, 8, ddG_data[4])
                    sheet.write(row, 9, ddG_data[5])
                    if calculate_delta_G_bind:
                        apo_delta_G_data = np.array(apo_lines[0][:-1].split(',')).astype(np.float64) / 2
                        apo_ddG_data = np.subtract(apo_delta_G_data, wt_apo_delta_G_dict.get(apo_state))
                        sheet.write(row, 2, float(ddG_data[0]) - float(apo_ddG_data[0]))
                        sheet.write(row, 4, float(ddG_data[1]) - float(apo_ddG_data[1]))
                        sheet.write(row, 7, float(ddG_data[3]) - float(apo_ddG_data[3]))
                    row += 1
                else:
                    wt_delta_G_data_1, wt_delta_G_data_2 = wt_delta_G_dict.get(nsp)
                    if calculate_delta_G_bind:
                        wt_apo_delta_G_data_1, wt_apo_delta_G_data_2 = wt_apo_delta_G_dict.get(apo_state)
                    sheet.write(row, 0, nsp + '_1')
                    delta_G_data_1 = np.array(lines[0][:-1].split(',')).astype(np.float64)
                    ddG_data_1 = np.subtract(delta_G_data_1, wt_delta_G_data_1)
                    sheet.write(row, 1, ddG_data_1[0])
                    sheet.write(row, 3, ddG_data_1[1])
                    sheet.write(row, 5, ddG_data_1[2])
                    sheet.write(row, 6, ddG_data_1[3])
                    sheet.write(row, 8, ddG_data_1[4])
                    sheet.write(row, 9, ddG_data_1[5])
                    if calculate_delta_G_bind:
                        apo_delta_G_data_1 = np.array(apo_lines[0][:-1].split(',')).astype(np.float64)
                        apo_ddG_data_1 = np.subtract(apo_delta_G_data_1, wt_apo_delta_G_data_1)
                        sheet.write(row, 2, float(ddG_data_1[0]) - float(apo_ddG_data_1[0]))
                        sheet.write(row, 4, float(ddG_data_1[1]) - float(apo_ddG_data_1[1]))
                        sheet.write(row, 7, float(ddG_data_1[3]) - float(apo_ddG_data_1[3]))
                    row += 1
                    sheet.write(row, 0, nsp + '_2')
                    delta_G_data_2 = np.array(lines[1][:-1].split(',')).astype(np.float64)
                    ddG_data_2 = np.subtract(delta_G_data_2, wt_delta_G_data_2)
                    sheet.write(row, 1, ddG_data_2[0])
                    sheet.write(row, 3, ddG_data_2[1])
                    sheet.write(row, 5, ddG_data_2[2])
                    sheet.write(row, 6, ddG_data_2[3])
                    sheet.write(row, 8, ddG_data_2[4])
                    sheet.write(row, 9, ddG_data_2[5])
                    if calculate_delta_G_bind:
                        apo_delta_G_data_2 = np.array(apo_lines[1][:-1].split(',')).astype(np.float64)
                        apo_ddG_data_2 = np.subtract(apo_delta_G_data_2, wt_apo_delta_G_data_2)
                        sheet.write(row, 2, float(ddG_data_2[0]) - float(apo_ddG_data_2[0]))
                        sheet.write(row, 4, float(ddG_data_2[1]) - float(apo_ddG_data_2[1]))
                        sheet.write(row, 7, float(ddG_data_2[3]) - float(apo_ddG_data_2[3]))
                    row += 1
    workbook.save(site + '.xls')

def calculate_inhibitor_ddG(inhibitor, site, native_aa):
    delta_G_data_file = inhibitor + '/' + site + '/' + site + '_'+ native_aa + '.dat'
    with open(delta_G_data_file) as pf:
        lines = pf.readlines()
    if len(lines) == 1:
        wt_inhibitor_delta_G_data = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
    else:
        wt_inhibitor_delta_G_data = (np.array(lines[0][:-1].split(',')).astype(np.float64), \
                np.array(lines[1][:-1].split(',')).astype(np.float64))
    apo_state = 'apo_' + inhibitor.split('_')[1]
    apo_delta_G_data_file = apo_state + '/' + site + '/' + site + '_'+ native_aa + '.dat'
    with open(apo_delta_G_data_file) as pf:
        lines = pf.readlines()
    if len(lines) == 1:
        wt_apo_delta_G_data = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
    else:
        wt_apo_delta_G_data = (np.array(lines[0][:-1].split(',')).astype(np.float64), \
                np.array(lines[1][:-1].split(',')).astype(np.float64))
    workbook_read = xlrd.open_workbook(inhibitor + '/' + site + '.xls')
    workbook = copy(workbook_read)
    for aa in ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', \
                'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']:
        mutation = site + aa
        # sheet_read = workbook_read.sheet_by_name(mutation)
        sheet = workbook.get_sheet(mutation)
        apo_delta_G_data_file = apo_state + '/' + site + '/' + site + '_'+ aa + '.dat'
        calculate_delta_G_bind = False
        if os.path.isfile(apo_delta_G_data_file):
            calculate_delta_G_bind = True
            with open(apo_delta_G_data_file) as pf:
                apo_lines = pf.readlines()
        delta_G_data_file = inhibitor + '/' + site + '/' + site + '_'+ aa + '.dat'
        if os.path.isfile(delta_G_data_file):
            with open(delta_G_data_file) as pf:
                lines = pf.readlines()
            if len(lines) == 1:
                inhibitor_delta_G_data = np.array(lines[0][:-1].split(',')).astype(np.float64) / 2
                inhibitor_ddG_data = np.subtract(inhibitor_delta_G_data, wt_inhibitor_delta_G_data)
                sheet.write(21, 0, inhibitor)
                sheet.write(21, 1, inhibitor_ddG_data[0])
                sheet.write(21, 3, inhibitor_ddG_data[1])
                sheet.write(21, 5, inhibitor_ddG_data[2])
                sheet.write(21, 6, inhibitor_ddG_data[3])
                sheet.write(21, 8, inhibitor_ddG_data[4])
                sheet.write(21, 9, inhibitor_ddG_data[5])
                if calculate_delta_G_bind:
                    apo_delta_G_data = np.array(apo_lines[0][:-1].split(',')).astype(np.float64) / 2
                    apo_ddG_data = np.subtract(apo_delta_G_data, wt_apo_delta_G_data)
                    sheet.write(21, 2, float(inhibitor_ddG_data[0]) - apo_ddG_data[0])
                    sheet.write(21, 4, float(inhibitor_ddG_data[1]) - apo_ddG_data[1])
                    sheet.write(21, 7, float(inhibitor_ddG_data[3]) - apo_ddG_data[3])
            else:
                inhibitor_delta_G_data_1 = np.array(lines[0][:-1].split(',')).astype(np.float64)
                inhibitor_ddG_data_1 = np.subtract(inhibitor_delta_G_data_1, wt_inhibitor_delta_G_data[0])
                sheet.write(21, 0, inhibitor + '_1')
                sheet.write(21, 1, inhibitor_ddG_data_1[0])
                sheet.write(21, 3, inhibitor_ddG_data_1[1])
                sheet.write(21, 5, inhibitor_ddG_data_1[2])
                sheet.write(21, 6, inhibitor_ddG_data_1[3])
                sheet.write(21, 8, inhibitor_ddG_data_1[4])
                sheet.write(21, 9, inhibitor_ddG_data_1[5])
                if calculate_delta_G_bind:
                    apo_delta_G_data_1 = np.array(apo_lines[0][:-1].split(',')).astype(np.float64)
                    apo_ddG_data_1 = np.subtract(apo_delta_G_data_1, wt_apo_delta_G_data[0])
                    sheet.write(21, 2, float(inhibitor_ddG_data_1[0]) - apo_ddG_data_1[0])
                    sheet.write(21, 4, float(inhibitor_ddG_data_1[1]) - apo_ddG_data_1[1])
                    sheet.write(21, 7, float(inhibitor_ddG_data_1[3]) - apo_ddG_data_1[3])
                inhibitor_delta_G_data_2 = np.array(lines[1][:-1].split(',')).astype(np.float64)
                inhibitor_ddG_data_2 = np.subtract(inhibitor_delta_G_data_2, wt_inhibitor_delta_G_data[1])
                sheet.write(21, 0, inhibitor + '_2')
                sheet.write(21, 1, inhibitor_ddG_data_2[0])
                sheet.write(21, 3, inhibitor_ddG_data_2[1])
                sheet.write(21, 5, inhibitor_ddG_data_2[2])
                sheet.write(21, 6, inhibitor_ddG_data_2[3])
                sheet.write(21, 8, inhibitor_ddG_data_2[4])
                sheet.write(21, 9, inhibitor_ddG_data_2[5])
                if calculate_delta_G_bind:
                    apo_delta_G_data_2 = np.array(apo_lines[1][:-1].split(',')).astype(np.float64)
                    apo_ddG_data_2 = np.subtract(apo_delta_G_data_2, wt_apo_delta_G_data[1])
                    sheet.write(21, 2, float(inhibitor_ddG_data_2[0]) - apo_ddG_data_2[0])
                    sheet.write(21, 4, float(inhibitor_ddG_data_2[1]) - apo_ddG_data_2[1])
                    sheet.write(21, 7, float(inhibitor_ddG_data_2[3]) - apo_ddG_data_2[3])
    workbook.save(inhibitor + '/' + site + '.xls')


if __name__ == '__main__':
    # delta_G_data_file: total[0],shell[1],substrate[2],residue[3],interaction[4],constraint[5]
    args = parse_arguments()
    site, native_aa = args.site_aa
    if not args.inhibitor:
        wt_apo_delta_G_data, wt_apo_delta_G_dict, wt_delta_G_dict = read_variant_delta_G(site, native_aa)
        calculate_variant_ddG(site, wt_apo_delta_G_data, wt_apo_delta_G_dict, wt_delta_G_dict)
    else:
        calculate_inhibitor_ddG(args.inhibitor, site, native_aa)
        
