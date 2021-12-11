
import numpy as np


#===========================================================================
# Print the summary of results at each iteration
#===========================================================================
def print_iteration_IPM(LO, coeff_matrix, rhs, iteration, norm_of_residual, is_sign_changed, alpha_star, alpha_hat, run_time):
	if iteration == 0:
		if LO.Params.LO_Verbosity > 0:
			print()

		if LO.Params.LO_Verbosity == 1 and LO.Params.Is_Quantum == True:
			print("{:18}{:}{:^16}{}".format(" ", "Objective", " ","Residual"), sep='', end="")

		elif LO.Params.LO_Verbosity == 1 and LO.Params.Is_Quantum  == False:
			print("{:18}{:}{:^16}{}".format(" ", "Objective", " ","Residual"), sep='', end="")

		elif LO.Params.LO_Verbosity == 2 and LO.Params.Is_Quantum  == True:
			print("{:18}{:}{:16}{:}{:24}{:^30s}".format(" ", "Objective", " ","Residual", " ", "Linear system" ), sep='', end="")
			if "II-" in LO.Params.Method: 
				print("{:2}{}".format("","Step-length" ), sep='', end="")

		elif LO.Params.LO_Verbosity == 2 and LO.Params.Is_Quantum  == False:
			print("{:18}{:}{:^16}{}".format(" ", "Objective", " ","Residual"), sep='', end="")
			if "II-" in LO.Params.Method:
				print("{:47}{}".format( " ", "Step-length" ), sep='', end="")


		if LO.Params.LO_Verbosity > 0:
			print()
			print("{:^4s}  ".format("Iter"), sep="", end="")
			print("{:^15s} ".format("Primal"), sep="", end="")
			print("{:^15s}  ".format("Dual"), sep="", end="")
			print("{:^8s} ".format("Primal"), sep="", end="")
			print("{:^8s}  ".format("Dual"), sep="", end="")
			print("{:^8s}".format("Compl"), sep="", end="")


		if LO.Params.LO_Verbosity > 0 and LO.Params.Is_Quantum  == True:
			print("  {:^8s}".format("||Resl||"), sep="", end="")


		# if LO.Params.LO_Verbosity == 2 and LO.Params.Is_Quantum  == True:
		# 	print(" {:^8s}".format("Prob"), sep="", end="")


		if LO.Params.LO_Verbosity == 2:
			print("  {:^8s} ".format("Cond-Num"), sep="", end="")
			print("{:^8s} ".format("||M||"), sep="", end="")
			print("{:^8s}  ".format("||RHS||"), sep="", end="")

			if "II-" in LO.Params.Method: 

				print("{:^8s} ".format("a-star"), sep="", end="")
				print("{:^8s}  ".format("a-hat"), sep="", end="")


		if LO.Params.LO_Verbosity > 0:
			print("{:>6s}".format("Time"), sep="", end="")
			print()
	


	if LO.Params.LO_Verbosity > 0:
		print("{:>4d}  ".format(iteration), sep="", end="")
		print("{:>15.8e} ".format(LO.primal_objective()), sep="", end="")
		print("{:>15.8e}  ".format(LO.dual_objective()), sep="", end="")
		print("{:>8.2e} ".format(LO.primal_residual()), sep="", end="")
		print("{:>8.2e} ".format(LO.dual_residual()), sep="", end="")
		if LO.complementarity() > 0:
			print(" {:>8.2e}".format(LO.complementarity()), sep="", end="")
		else:
			print(" {:>7.2e}".format(LO.complementarity()), sep="", end="")



	
	if LO.Params.LO_Verbosity > 0 and LO.Params.Is_Quantum  == True:
		if is_sign_changed == True:
			print(" *{:>8.2e}".format(norm_of_residual), sep="", end="")
		else:
			print("  {:>8.2e}".format(norm_of_residual), sep="", end="")

	
	
	# if LO.Params.LO_Verbosity == 2 and LO.Params.Is_Quantum  == True:
	# 	print(" {:>8.2e}".format(probability), sep="", end="")

	
	
	if LO.Params.LO_Verbosity == 2:
		condition_number 		= np.linalg.cond(coeff_matrix)
		print("  {:>8.2e} ".format(condition_number), sep="", end="")
		print("{:>8.2e} ".format(np.linalg.norm(coeff_matrix, 2)), sep="", end="")
		print("{:>8.2e}  ".format(np.linalg.norm(rhs) ), sep="", end="")

		if "II-" in LO.Params.Method: 

			print("{:>8.2e} ".format(alpha_star), sep="", end="")
			print("{:>8.2e}  ".format(alpha_hat), sep="", end="")

	
	if LO.Params.LO_Verbosity > 0:
		print("{:>5.0f}s".format(run_time), sep="", end="")
		print()



#===========================================================================
# Print the final results
#===========================================================================
def print_final_IPM(LO, iteration, run_time):
	precision_digit = - int(np.log10(LO.Params.LO_Precision)) + 2

	if LO.Params.LO_Verbosity > -1:
		print()
		print("The algorithm stopped after {:d} iterations in {:.2f} seconds.".format(iteration, run_time))


	if LO.Params.LO_Verbosity == 0:
		print()
		print("{:20}{:<15.8e}".format("Primal objective:", LO.primal_objective()))
		print("{:20}{:<15.8e}".format("Dual objective:", LO.dual_objective()))
		print()
		print("{:20}{:<8.2e}".format("Primal residual:", LO.primal_residual()))
		print("{:20}{:<8.2e}".format("Dual residual:", LO.dual_residual()))
		print("{:20}{:<8.2e}".format("Complementraty:", LO.complementarity()))


	if LO.Params.LO_Verbosity > 0:
		np.set_printoptions(precision=precision_digit, formatter={'float': ('{:+.'+str(precision_digit)+'f}').format})
		print()
		print("{:20}".format("Primal variables:"), LO.x, sep="")
		print("{:20}".format("Dual slacks:"), LO.s, sep="")
		print("{:20}".format("Dual variables:"), LO.y, sep="")
		print()
		print("{:20}{:<15.8e}".format("Primal objective:", LO.primal_objective()))
		print("{:20}{:<15.8e}".format("Dual objective:", LO.dual_objective()))
		print()
		print("{:20}{:<8.2e}".format("Primal residual:", LO.primal_residual()))
		print("{:20}{:<8.2e}".format("Dual residual:", LO.dual_residual()))
		print("{:20}{:<8.2e}".format("Complementraty:", LO.complementarity()))
		print()
		print("{:20}".format("Number of Iter:"), iteration, sep="")
		print("{:20}{:<.2f}".format("Run time:", run_time), sep="")






