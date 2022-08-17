# install forgi:
# pip install forgi

import csv
import forgi
import forgi.graph.bulge_graph as fgb


with open('ev_DesignMatrix_SeqPlusDotBracket.csv', 'r') as csv_input:
	with open('ev_DesignMatrix_SeqPlusDotBracket_forgi.csv', 'w') as csv_output:
		readCSV = csv.reader(csv_input, delimiter=',')
		writeCSV = csv.writer(csv_output, lineterminator='\n', delimiter=',')
		# all_col = []
		row = next(readCSV)
		seq_index = row.index('seq')
		dotbracket_index = row.index('dotbracket')
		row.append('element_string')
		row.append('element_string_number')
		writeCSV.writerow(row)
		# all_col.append(row)
		for row in readCSV:
			# print(row[0])
			seq = row[seq_index]
			# print('seq')
			# print(seq)
			dotbracket = row[dotbracket_index]
			# print("dotbracket:")
			# print(dotbracket)
			if len(seq) != len(dotbracket):
				print('Different Lengths, failed:')
				print(row)
				break
			tmp_file = open("tmp.txt","w")
			tmp_file.write(seq + '\n')
			tmp_file.write(dotbracket + '\n')
			tmp_file.close()
			cg = forgi.load_rna('tmp.txt', allow_many=False)
			result = fgb.BulgeGraph.to_element_string(cg, with_numbers=True)
			result = result.splitlines()
			row.append(str(result[0]))
			row.append(str(result[1]))
			writeCSV.writerow(row)
			# all_col.append(row)
		# writeCSV.writerows(all_col)



# bg = fgb.BulgeGraph.from_dotbracket('((..))..((..))')

# rna_seq_dotb = ['CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUU', '((((((((((..(((((((.......)))))))......).(((.((.......))))))..)))))))).']
# tmp_file = open("tmp.txt","w")
# tmp_file.write(rna_seq_dotb[0]+ '\n')
# tmp_file.write(rna_seq_dotb[1]+ '\n')
# tmp_file.writelines(rna_seq_dotb)
# tmp_file.close()
# cg = forgi.load_rna('tmp.txt', allow_many=False)
# cg = forgi.load_rna("example/testdotbracket", allow_many=False)
# fgb.BulgeGraph.to_dotbracket_string(cg)
# result = fgb.BulgeGraph.to_element_string(cg, with_numbers=True)
# result = result.splitlines()

# print("len:")
# print(len(result))
# print(result[0])
# print(result[1])

