List of files:
- Force.m
- LHS.m
- main_gen.m
- mesh.m
- monomials.m
- plot.m
- RHS.m
- RHS_fromData.m
- Sol.m
- sourcePnts.m

Description:
- Force.m 
  Calculates the source function on each time step for given point
- LHS.m
  Assembles the matrix for linear system
- main_gen.m
  Main file 
- mesh.m
  Generates domain/boundary/source points in case if data is not given and one's testing true solution
- monomials.m
  Calculates values of monomials (up to first 10)
- plot_res.m
  Visalization of results
- RHS.m
  Assembles the right-hand side vector at each time step in case of analytical source function
- RHS_fromData.m
  Assembles the right-hand side vector at each time step in case of given data
- Sol.m
  Calculates the solution (temperature) at each time step for given point
- sourcePnts.m
  Generates source points in case of given data
  	
To get the results, run main_gen.m file and then plot_res.m.
