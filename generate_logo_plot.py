import argparse
import logomaker
import numpy as np
import pandas as pd
import statistics
import xlrd

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('xls', type=str)
    parser.add_argument('-asymm', action="store_true")
    parser.add_argument('-type', '--fitness_type', type=str, choices=['tolerable', 'escape'], required=True)
    return parser.parse_args()

def calculate_tolerable_fitness(row_values, fitness_type):
    mutation = row_values[0]
    ddG_fold = row_values[1]
    ddG_bind_inhibitor = row_values[2]
    ddG_sub_inhibitor = row_values[3]
    ddG_bind_inhibitor_zscore = row_values[4]
    ddG_bind_sub = np.array(row_values[5:]).astype(np.float64).reshape(-1, 2)
    ddG_bind = ddG_bind_sub[:, 0]
    ddG_sub = ddG_bind_sub[:, 1]
    # max_ddG_bind, second_max_ddG_bind = sorted(ddG_bind, reverse=True)[:2]
    average_ddG_bind = statistics.mean(ddG_bind)
    # max_ddG_sub, second_max_ddG_sub = sorted(ddG_sub, reverse=True)[:2]
    average_ddG_sub = statistics.mean(ddG_sub)
    if fitness_type == 'tolerable':
        # fitness = pow(1.2, -(ddG_fold + average_ddG_bind + (max_ddG_bind + second_max_ddG_bind)/2 \
        #         + average_ddG_sub + (max_ddG_sub + second_max_ddG_sub)/2))
        fitness = pow(1.2, -(ddG_fold + 2 * average_ddG_bind + 2 * average_ddG_sub))
    elif fitness_type == 'escape':
        # fitness = pow(1.2, -(ddG_fold \
        #         + average_ddG_bind + (max_ddG_bind + second_max_ddG_bind)/2 - ddG_bind_inhibitor - ddG_bind_inhibitor_zscore \
        #         + average_ddG_sub + (max_ddG_sub + second_max_ddG_sub)/2 - ddG_sub_inhibitor))
        fitness = pow(1.2, -(ddG_fold + average_ddG_bind + average_ddG_sub \
                + average_ddG_bind - ddG_bind_inhibitor \
                + average_ddG_bind - ddG_bind_inhibitor_zscore \
                + average_ddG_sub - ddG_sub_inhibitor))
    return mutation, fitness

def convert_xls_to_csv(workbook_name, chains, fitness_type):
    workbook_read = xlrd.open_workbook(workbook_name)
    for chain in range(chains):
        sheet_read = workbook_read.sheet_by_index(chain)
        positions = np.empty(0)
        with open(workbook_name[:-4] + '.' + fitness_type + '.' + str(chain + 1) + 'pssm', 'w') as pf:
            pf.write('pos	A	G	I	L	P	V	F	W	Y	D	E	R	H	K	S	T	C	M	N	Q\n1	')
            current_pos = sheet_read.cell_value(1, 0)[:-1]
            npos = 1
            for row in range(1, sheet_read.nrows):
                row_values = sheet_read.row_values(row)
                mutation, probability = calculate_tolerable_fitness(row_values, fitness_type)
                if mutation[:-1] == current_pos:
                    pf.write(str(round(probability, 3)) + '	')
                else:
                    positions = np.append(positions, current_pos)
                    current_pos = mutation[:-1]
                    npos += 1
                    pf.write('\n' + str(npos) + '	' + str(round(probability, 3)) + '	')
            positions = np.append(positions, current_pos)
            pf.write('\n')
    return positions

def logo_plot(workbook_name, chains, positions, fitness_type):
    for chain in range(chains):
        df = pd.read_csv(workbook_name[:-4] + '.' + fitness_type + '.' + str(chain + 1) + 'pssm', delim_whitespace=True, index_col=0)
        logo = logomaker.Logo(df, font_name='Arial Rounded MT Bold', figsize=(len(positions), 10))
        logo.ax.set_xticks(range(1, len(positions) + 1))
        logo.ax.set_xticklabels(positions)
        # logo.fig.show()
        logo.fig.savefig(workbook_name[:-4] + '.' + fitness_type + '.' + str(chain + 1) + 'png', bbox_inches='tight', dpi=200)

if __name__ == '__main__':
    args = parse_arguments()
    chains = 2 if args.asymm else 1
    positions = convert_xls_to_csv(args.xls, chains, args.fitness_type)
    logo_plot(args.xls, chains, positions, args.fitness_type)
