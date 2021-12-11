import xlwt 
import os.path
from os import walk
import os 

files 		= [filename for filename in sorted(os.listdir(os.getcwd())) if (not os.path.isdir(filename) and filename.split("_")[0] == 'runner')]


book 		= xlwt.Workbook(encoding="utf-8")
sheet1 		= book.add_sheet("detailed")
sheet2		= book.add_sheet("summary")

headers 	= ["File name", "Seed", "Norm of RHS vector", "Norm of matrix", "Condition number", "Num of ancillae qubits", "Expansion order", "Num of time slices" , "Norm of difference", "Norm of residual", "Is sign changed", "Probability", "Time (s)", "Circuit depth", "Circuit width", "Norm of RHS M-vector", "Norm of M-matrix", "M-condition number", "M-norm of difference", "M-norm of residual"]
summary 	= ["File name", "Seed",  "Norm of RHS vector", "Norm of matrix", "Condition number", "Norm of difference", "Norm of residual", "Is sign changed", "Probability", "Time (s)", "Circuit depth", "Circuit width"]


for (col, header) in enumerate(headers):	
	sheet1.write(0, col, header)


for (col, header) in enumerate(summary):	
	sheet2.write(0, col, header)



def getValue(lines, header, file_name):

	if header == "File name":
		result 	= file_name

	elif header == "Seed":
		result 	= int(filename.split("_")[3][4:])

	else: 
		for line in reversed(lines):
			if line.split(":")[0].strip() == header:
				try: 	result 	= float(line.split(":")[1].strip())
				except: result 	= line.split(":")[1].strip()
				break

	return result



for (row, filename) in enumerate(files):
	with open(filename) as f:
		
		lines 	= f.read().splitlines()

		if lines[-1].strip() == 61*"=":

			for (col, header) in enumerate(headers):
				val 		= getValue(lines, header, filename)
				sheet1.write(row+1, col, val)
			
			for (col, header) in enumerate(summary):
				val 		= getValue(lines, header, filename)
				sheet2.write(row+1, col, val)
				 
			print ("Read file:",filename)

book.save("results.xls")	


