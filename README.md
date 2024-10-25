There are detailed comments related to the following:

Init.c -> MGInitialize()
ibm_io.c-> ibm_read_Icem()
ibm.c-> ibm_search_advanced()
ibm.c->ibm_search_advanced_rev()
ibm.c->point_cell_advanced*()
bcs.c->InflowFlux()
bcs.c->OutflowFlux()
bcs.c->FormBCS()


*multiple variations of this functions are present in the code.
There are detailed comments related to the following:

Init.c -> MGInitialize()
ibm_io.c-> ibm_read_Icem()
ibm.c-> ibm_search_advanced()
ibm.c->ibm_search_advanced_rev()
ibm.c->point_cell_advanced*()
bcs.c->InflowFlux()
bcs.c->OutflowFlux()
bcs.c->FormBCS()

*multiple variations of this functions are present in the code.

-------------------------------------------------------------------------
-------------------------------------------------------------------------
Flag notes: 

/* immeresed value>1 determines the type of correction
   1      constant velocity correction
   2      proportional (needs FluxOutSum>0)
   3      proportional to flux
   4      proportional to normal velocity (flux/area)
*/
