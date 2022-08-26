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

def read_delta_G_from_file(delta_G_data_file):
    if os.path.isfile(delta_G_data_file):
        with open(delta_G_data_file) as pf:
            lines = pf.readlines()
        if len(lines) == 1:
            delta_G_data = np.array([lines[0][:-1].split(',')]).astype(np.float64) / 2
        else:
            delta_G_data = np.array([lines[0][:-1].split(','), lines[1][:-1].split(',')]).astype(np.float64)
    else:
        delta_G_data = np.zeros((0, 6))
    return delta_G_data

def read_variant_delta_G(site, native_aa):
    # Read wild type protease-nsp complex delta G and wild type apo protease delta G
    wt_delta_G_dict = dict()
    for nsp in ['nsp4-nsp5_7T70', 'nsp5-nsp6_7T8M', 'nsp6-nsp7_7MB6', 'nsp7-nsp8_7T8R', \
            'nsp8-nsp9_7T9Y', 'nsp9-nsp10_7TA4', 'nsp10-nsp11_7TA7', 'nsp12-nsp13_7TB2', \
            'nsp13-nsp14_7TBT', 'nsp14-nsp15_7TA4', 'nsp15-nsp16_7TC4', 'apo_6YB7']:
        wt_delta_G_data = read_delta_G_from_file(nsp + '/' + site + '/' + site + '_'+ native_aa + '.dat')
        wt_delta_G_dict[nsp] = wt_delta_G_data
        if not nsp.startswith('apo_'):
            apo_state = 'apo_' + nsp.split('_')[1]
            wt_apo_delta_G_data = read_delta_G_from_file(apo_state + '/' + site + '/' + site + '_'+ native_aa + '.dat')
            wt_delta_G_dict[apo_state] = wt_apo_delta_G_data
    return wt_delta_G_dict

def calculate_variant_ddG(site, wt_delta_G_dict):
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
        # Declare ddG fold data array
        ddG_fold = np.zeros(6)
        ddG_fold_n = 0
        row = 2
        # Calculate protease ddG for each amino acid type
        for nsp in ['nsp4-nsp5_7T70', 'nsp5-nsp6_7T8M', 'nsp6-nsp7_7MB6', 'nsp7-nsp8_7T8R', \
                'nsp8-nsp9_7T9Y', 'nsp9-nsp10_7TA4', 'nsp10-nsp11_7TA7', 'nsp12-nsp13_7TB2', \
                'nsp13-nsp14_7TBT', 'nsp14-nsp15_7TA4', 'nsp15-nsp16_7TC4', 'apo_6YB7']:
            # Calculate holo protease ddG total for each amino acid type
            if not nsp.startswith('apo_'):
                delta_G_data = read_delta_G_from_file(nsp + '/' + site + '/' + site + '_'+ aa + '.dat')
                wt_delta_G_data = wt_delta_G_dict.get(nsp)
                ddG_data = np.zeros((0, 6))
                for chain in range(min(delta_G_data.shape[0], wt_delta_G_data.shape[0])):
                    ddG = np.subtract(delta_G_data[chain], wt_delta_G_data[chain])
                    ddG_data = np.append(ddG_data, [ddG], axis=0)
                    sheet.write(row + chain, 0, nsp)
                    sheet.write(row + chain, 1, ddG[0])
                    sheet.write(row + chain, 3, ddG[1])
                    sheet.write(row + chain, 5, ddG[2])
                    sheet.write(row + chain, 6, ddG[3])
                    sheet.write(row + chain, 8, ddG[4])
                    sheet.write(row + chain, 9, ddG[5])
            # Calculate apo protease ddG fold and ddG bind for each amino acid type
            apo_state = 'apo_' + nsp.split('_')[1]
            apo_delta_G_data = read_delta_G_from_file(apo_state + '/' + site + '/' + site + '_'+ aa + '.dat')
            wt_apo_delta_G_data = wt_delta_G_dict.get(apo_state)
            for chain in range(min(apo_delta_G_data.shape[0], wt_apo_delta_G_data.shape[0])):
                apo_ddG = np.subtract(apo_delta_G_data[chain], wt_apo_delta_G_data[chain])
                ddG_fold = np.add(ddG_fold, apo_ddG)
                ddG_fold_n += 1
                if ddG_data.shape[0] > 0 and not nsp.startswith('apo_'):
                        sheet.write(row + chain, 2, float(ddG_data[chain][0]) - float(apo_ddG[0]))
                        sheet.write(row + chain, 4, float(ddG_data[chain][1]) - float(apo_ddG[1]))
                        sheet.write(row + chain, 7, float(ddG_data[chain][3]) - float(apo_ddG[3]))
            if not nsp.startswith('apo_'):
                row += ddG_data.shape[0]
        # Write averaged apo protease ddG
        if ddG_fold_n > 0:
            ddG_fold = ddG_fold / ddG_fold_n
            sheet.write(1, 0, 'apo')
            sheet.write(1, 1, ddG_fold[0])
            sheet.write(1, 2, 'NA')
            sheet.write(1, 3, ddG_fold[1])
            sheet.write(1, 4, 'NA')
            sheet.write(1, 5, ddG_fold[2])
            sheet.write(1, 6, ddG_fold[3])
            sheet.write(1, 7, 'NA')
            sheet.write(1, 8, ddG_fold[4])
            sheet.write(1, 9, ddG_fold[5])
    workbook.save(site + '.xls')

def calculate_inhibitor_ddG(inhibitor, site, native_aa):
    wt_inhibitor_delta_G_data = read_delta_G_from_file(inhibitor + '/' + site + '/' + site + '_'+ native_aa + '.dat')
    apo_state = 'apo_' + inhibitor.split('_')[1]
    wt_apo_delta_G_data = read_delta_G_from_file(apo_state + '/' + site + '/' + site + '_'+ native_aa + '.dat')
    workbook_read = xlrd.open_workbook(site + '.xls')
    workbook = copy(workbook_read)
    for aa in ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', \
                'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']:
        mutation = site + aa
        sheet = workbook.get_sheet(mutation)
        apo_ddG_data = np.zeros((0, 6))
        apo_delta_G_data = read_delta_G_from_file(apo_state + '/' + site + '/' + site + '_'+ aa + '.dat')
        for chain in range(min(apo_delta_G_data.shape[0], wt_apo_delta_G_data.shape[0])):
            apo_ddG = np.subtract(apo_delta_G_data[chain], wt_apo_delta_G_data[chain])
            apo_ddG_data = np.append(apo_ddG_data, [apo_ddG], axis=0)
        inhibitor_delta_G_data = read_delta_G_from_file(inhibitor + '/' + site + '/' + site + '_'+ aa + '.dat')
        for chain in range(min(inhibitor_delta_G_data.shape[0], wt_inhibitor_delta_G_data.shape[0])):
            inhibitor_ddG = np.subtract(inhibitor_delta_G_data[chain], wt_inhibitor_delta_G_data[chain])
            sheet.write(21 + chain, 0, inhibitor)
            sheet.write(21 + chain, 1, inhibitor_ddG[0])
            sheet.write(21 + chain, 3, inhibitor_ddG[1])
            sheet.write(21 + chain, 5, inhibitor_ddG[2])
            sheet.write(21 + chain, 6, inhibitor_ddG[3])
            sheet.write(21 + chain, 8, inhibitor_ddG[4])
            sheet.write(21 + chain, 9, inhibitor_ddG[5])
            if apo_ddG_data.shape[0] > 0:
                sheet.write(21 + chain, 2, float(inhibitor_ddG[0]) - apo_ddG_data[chain][0])
                sheet.write(21 + chain, 4, float(inhibitor_ddG[1]) - apo_ddG_data[chain][1])
                sheet.write(21 + chain, 7, float(inhibitor_ddG[3]) - apo_ddG_data[chain][3])
    workbook.save(inhibitor + '/' + site + '.xls')


if __name__ == '__main__':
    # delta_G_data_file: total[0],shell[1],substrate[2],residue[3],interaction[4],constraint[5]
    args = parse_arguments()
    site, native_aa = args.site_aa
    if not args.inhibitor:
        wt_delta_G_dict = read_variant_delta_G(site, native_aa)
        calculate_variant_ddG(site, wt_delta_G_dict)
    else:
        calculate_inhibitor_ddG(args.inhibitor, site, native_aa)
