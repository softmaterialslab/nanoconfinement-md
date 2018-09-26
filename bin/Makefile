# This submit a job in bigred2

JOBSCR = iu_cluster_job_script.pbs
TEST = test.pbs

submit: 
	@echo "Launching the job";
	qsub $(JOBSCR)

test: 
	@echo "Launching the test job for nanoconfinement-md";
	qsub $(TEST)