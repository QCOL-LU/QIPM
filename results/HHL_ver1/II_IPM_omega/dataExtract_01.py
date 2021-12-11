import xlwt 
import os.path
from os import walk
import os 

files 		= [filename for filename in sorted(os.listdir(os.getcwd())) if (not os.path.isdir(filename) and filename.split("_")[0] == 'runner')]


book 		= xlwt.Workbook(encoding="utf-8")
sheet1 		= book.add_sheet("detailed")
sheet2		= book.add_sheet("summary")

headers 	= ["File name", "Seed", "Omega", "Number of variables", "Number of constraints", "Norm of primal sol", "Norm of dual sol", "Norm of dual slacks" , "Norm of vector b", "Norm of vector c", "Norm of matrix A", "Condition number", "Optimal objective", "Primal objective", "Dual objective", "Primal residual", "Dual residual", "Complementraty", "Number of Iter", "Run time"]
summary 	= ["File name", "Seed", "Omega", "Number of variables", "Number of constraints", "Norm of vector b", "Norm of matrix A", "Condition number", "Optimal objective", "Primal objective", "Dual objective", "Primal residual", "Dual residual", "Complementraty", "Number of Iter", "Run time"]

for (col, header) in enumerate(headers):	
	sheet1.write(0, col, header)


for (col, header) in enumerate(summary):	
	sheet2.write(0, col, header)



def getValue(lines, header, file_name):
	result = False
	if header == "File name":
		result 	= file_name

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
		
		flag = 0
		for (col, header) in enumerate(headers):
			val 		= getValue(lines, header, filename)
			sheet1.write(row+1, col, val)
			
			if val == -1 and flag == 0: 
				flag = 1
				print(filename)
		
		for (col, header) in enumerate(summary):
			val 		= getValue(lines, header, filename)
			sheet2.write(row+1, col, val)
			 
		# print ("Read file:",filename)

book.save("results.xls")	


