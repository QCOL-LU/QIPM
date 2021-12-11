import numpy as np


#===========================================================================
# Print the summary of results at each iteration
#===========================================================================
def print_iteration_iterative_refinement(LO, r, r_bar, nabla, iteration, run_time):
	if iteration == 0 or LO.Params.LO_Verbosity > -1:
		if LO.Params.IR_Verbosity > 0:
			print()
		

		if LO.Params.IR_Verbosity == 1:
			if LO.Params.LO_Verbosity > -1: print(75*"-")
			print("{:18}{:}{:^16}{}".format(" ", "Objective", " ","Residual"), sep='', end="")
			

		elif LO.Params.IR_Verbosity == 2:
			if LO.Params.LO_Verbosity > -1: print(103*"-")
			print("{:18}{:}{:^16}{}{:17}{:^25s}".format(" ", "Objective", " ","Residual", " ", "Step-length" ), sep='', end="")
			


		if LO.Params.IR_Verbosity > 0:
			print()
			print("{:^4s}  ".format("Iter"), sep="", end="")
			print("{:^15s} ".format("Primal"), sep="", end="")
			print("{:^15s}  ".format("Dual"), sep="", end="")
			print("{:^8s} ".format("Primal"), sep="", end="")
			print("{:^8s}  ".format("Dual"), sep="", end="")
			print("{:^8s}".format("Compl"), sep="", end="")


		if LO.Params.IR_Verbosity == 2:
			print("  {:^8s} ".format("r error"), sep="", end="")
			print("{:^8s} ".format("r-bar"), sep="", end="")
			print("{:^8s}".format("nabla"), sep="", end="")


		if LO.Params.IR_Verbosity > 0:
			print(" {:>6s}".format("Time"), sep="", end="")
			print()
	


	if LO.Params.IR_Verbosity > 0:
		print("{:>4d}  ".format(iteration), sep="", end="")
		print("{:>15.8e} ".format(LO.primal_objective()), sep="", end="")
		print("{:>15.8e}  ".format(LO.dual_objective()), sep="", end="")
		print("{:>8.2e} ".format(LO.primal_residual()), sep="", end="")
		print("{:>8.2e}  ".format(LO.dual_residual()), sep="", end="")
		print("{:>8.2e}".format(LO.complementarity()), sep="", end="")


	
	if LO.Params.IR_Verbosity == 2:
		print("  {:>8.2e} ".format(r), sep="", end="")

		print("{:>8.2e} ".format(r_bar), sep="", end="")
		print("{:>8.2e}".format(nabla), sep="", end="")

	
	if LO.Params.IR_Verbosity > 0:
		print(" {:>5.0f}s".format(run_time), sep="", end="")
		print()



#===========================================================================
# Print the final results
#===========================================================================
def print_final_iterative_refinement(LO, iteration, run_time):
	precision_digit = - int(np.log10(LO.Params.IR_Precision)) + 2

	if LO.Params.IR_Verbosity > -1:
		print()
		print("The Iterative Refinement algorithm stopped after {:d} iterations in {:.2f} seconds.".format(iteration+1, run_time))


	if LO.Params.IR_Verbosity == 0:
		print()
		print("{:20}{:<15.8e}".format("Primal objective:", LO.primal_objective()))
		print("{:20}{:<15.8e}".format("Dual objective:", LO.dual_objective()))
		print()
		print("{:20}{:<8.2e}".format("Primal residual:", LO.primal_residual()))
		print("{:20}{:<8.2e}".format("Dual residual:", LO.dual_residual()))
		print("{:20}{:<8.2e}".format("Complementraty:", LO.complementarity()))


	if LO.Params.IR_Verbosity > 0:
		np.set_printoptions(precision=precision_digit, formatter={'float': ('{:+.'+str(precision_digit)+'f}').format})
		print()
		print("{:20}{:<15.8e}".format("Primal objective:", LO.primal_objective()))
		print("{:20}{:<15.8e}".format("Dual objective:", LO.dual_objective()))
		print()
		print("{:20}{:<8.2e}".format("Primal residual:", LO.primal_residual()))
		print("{:20}{:<8.2e}".format("Dual residual:", LO.dual_residual()))
		print("{:20}{:<8.2e}".format("Complementraty:", LO.complementarity()))
		print()
		print("{:20}".format("Primal variables:"), LO.x, sep="")
		print("{:20}".format("Dual slacks:"), LO.s, sep="")
		print("{:20}".format("Dual variables:"), LO.y, sep="")
		print()
		print("{:20}".format("Number of Iter:"), iteration, sep="")
		print("{:20}{:<.2f}".format("Run time:", run_time), sep="")





