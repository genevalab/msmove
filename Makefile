# Builds all the projects in the solution...
.PHONY: all_projects
all_projects: msmove 

# Builds project 'msmove'...
.PHONY: msmove
msmove: 
	make --directory="." --file=msmove.makefile

# Cleans all projects...
.PHONY: clean
clean:
	make --directory="." --file=msmove.makefile clean

