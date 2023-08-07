"# TwoFlankABModel" 
"# TwoFlankABModel" 
This package was developed by Dominic T Robson and Andreas C W Baas for the simulation of barchan swarms.

Contained within this source code are five python modules and an example of a runfile you might use to perform simulations using the model.  Other prerequisite packages you will require are numpy, matplotlib, mesa, and shapely (use version 1.8.0 as updated version currently contains some bugs)

The five python modules included in this source code are:
	ABModel.py - contains definitions for the Swarm (model) and Barchan (agent)
			and key functions relating to how the agents and model behave
	CollisionStuff.py - contains functions used to determine which agents collide and
				the outcome of those collisions
	FluxAndAvalanching.py - contains functions used for calculating the flux absorbed by 				a dune as well as definitions of the geometry of the objects
	GammaStuff.py - contains a function for calculating the maximum asymmetry threshold
	WinddirectionStuff.py - contains functions for converting between wind angles and 
				and vectors as well as for generating stochastic or 					predetermined wind regimes

Together these five contain everything you need to execute the Two-Flank Agent-Based model (TFABM), designed for simulating swarms of barchan dunes.

We have also provided an example of a runfile which walks the user through the basics of using the model and also provides an example of a simulation of a small swarm containing three dunes.







