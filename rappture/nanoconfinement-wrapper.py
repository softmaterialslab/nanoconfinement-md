# ----------------------------------------------------------------------
#  Python wrapper for Nanoconfinement code.
#
#  This wrapper script calls Rappture API to fetch input values and
#   passes them to nanoconfinement tool. The wrapper can be extended
#   to do higher level data manipulation to marshall user inputs.
#
# ----------------------------------------------------------------------

import Rappture
import sys, os, commands, string, shutil, math

# open the XML file containing the run parameters
driver = Rappture.library(sys.argv[1])

# Parse the rappture generated XML file to extract user input values
io = Rappture.PyXml(sys.argv[1])

salt_concentration = io['input.group(physical).number(salt_concentration).current'].value
print "salt concentration is %s" % salt_concentration

positive_valency = io['input.group(physical).integer(positive_valency).current'].value
print "positive valency is %s" % positive_valency

negative_valency = io['input.group(physical).integer(negative_valency).current'].value
print "negative_valency  s %s" % negative_valency

confinement_length = io['input.group(physical).number(confinement_length).current'].value
print "confinement_length is %s" % confinement_length

ion_diameter = io['input.group(physical).number(ion_diameter).current'].value
print "ion_diameter is %s" % ion_diameter

simulation_steps = io['input.group(computing).integer(simulation_steps).current'].value
print "simulation_steps is %s" % simulation_steps

simulation_params="_%.2f" % float(confinement_length)+"_%d" % int(positive_valency)+"_%d" % int(negative_valency)+"_%.2f" % float(salt_concentration)+"_%.3f" % float(ion_diameter)+"_%d" % int(simulation_steps);

shutil.rmtree('data',True)
if not os.path.exists('data'):
    os.makedirs('data')

os.system("use boost-1.62.0-mpich2-1.3-gnu-4.7.2")

runName='nanoconfine'

mpi_processors=round((int(simulation_steps) + 333333)/333333)
total_processors=mpi_processors*16

walltime = round(16*(1+(5.25 * math.exp(-mpi_processors/1.78))))

print "walltime is %d" % walltime
print "mpi_processors is %d" % mpi_processors
print "total_processors is %d" % total_processors

try:
     #exitStatus,stdOutput,stdError = Rappture.tools.executeCommand(
     #   ['mpirun','-np','1','md_simulation_confined_ions', '-Z', confinement_length, '-p', positive_valency, '-n', negative_valency, '-c',
     #   salt_concentration, '-d', ion_diameter, '-S', simulation_steps, '-f', simulation_params, '-v', 'false'], streamOutput=True)
	 
	 exitStatus,stdOutput,stdError = Rappture.tools.executeCommand(
	 ['submit','--venue','standby@conte','-w',walltime,'-n',total_processors, '-N','16', '--runName',runName, '--tailStdout', '--inputfile','data', 'nanoconfinement-r21',
		 '-Z', confinement_length, '-p', positive_valency, '-n', negative_valency, '-c', salt_concentration, 
		 '-d', ion_diameter, '-S', simulation_steps, '-f', simulation_params, '-v', 'false'], streamOutput=True)
 		 
except:
    sys.stderr.write('Error during execution of md_simulation_confined_ions')
    sys.exit(1);
	
# Reading standard output from the file
try:
	fid = open(runName+'.stdout','r')
	info = fid.readlines()
	fid.close()
	os.remove(runName+'.stdout') 
except:
	sys.stderr.write('Can not find the .stdout file')
	sys.exit(1);
# Setting standard output to GUI	
io['output.log']='\n'.join(info)
		
# Label output graph with title, x-axis label,
# Positive density profile
io['output.curve(positive_ion_density).about.label']='Density of Positive Ions'
io['output.curve(positive_ion_density).about.description']='Distribution of positive ions confined within the nanoparticle surfaces'
io['output.curve(positive_ion_density).xaxis.label']='z'
io['output.curve(positive_ion_density).xaxis.description']='Distance measure between the two Surfaces (z = 0 is the midpoint)'
io['output.curve(positive_ion_density).xaxis.units']='nm'
io['output.curve(positive_ion_density).yaxis.label']='Density'
io['output.curve(positive_ion_density).yaxis.description']='Density distribution of ions'
io['output.curve(positive_ion_density).yaxis.units']='M'

try:
	fid = open('data/p_density_profile'+simulation_params+'.dat','r')
	info = fid.readlines()
	fid.close()
	os.remove('data/p_density_profile'+simulation_params+'.dat') 
except:
	sys.stderr.write('Can not find the positive density results file')
	sys.exit(1);

# add density profile to xy data
xList = []
yList = []
for line in info:
	proLine=" ".join(line.split())
	d,m,e = proLine.split()
	xList.append(float(d))
	yList.append(float(m))

io['output.curve(positive_ion_density).component.xy']=(xList, yList)

# Negative density profile
io['output.curve(negative_ion_density).about.label']='Density of Negative Ions'
io['output.curve(negative_ion_density).about.description']='Distribution of negative ions confined within the nanoparticle surfaces'
io['output.curve(negative_ion_density).xaxis.label']='z'
io['output.curve(negative_ion_density).xaxis.description']='Distance measure between the two Surfaces (z = 0 is the midpoint)'
io['output.curve(negative_ion_density).xaxis.units']='nm'
io['output.curve(negative_ion_density).yaxis.label']='Density'
io['output.curve(negative_ion_density).yaxis.description']='Density distribution of ions'
io['output.curve(negative_ion_density).yaxis.units']='M'

try:
	fid = open('data/n_density_profile'+simulation_params+'.dat','r')
	info = fid.readlines()
	fid.close()
	os.remove('data/n_density_profile'+simulation_params+'.dat') 
except:
	sys.stderr.write('Can not find the negative density results file')
	sys.exit(1);
	
#remove the data folder
#os.rmdir('data')
shutil.rmtree('data')
		
# add density profile to xy data
xListN = []
yListN = []
for line in info:
	proLine=" ".join(line.split())
	d,m,e = proLine.split()
	xListN.append(float(d))
	yListN.append(float(m))
	
io['output.curve(negative_ion_density).component.xy']=(xListN, yListN)

#multi curve for combined +/- density profiles
io['output.curve(positive_ion_density_Combined).about.label']='Density of Positive Ions'
io['output.curve(positive_ion_density_Combined).about.group']='Density of Positive/Negative Ions'
io['output.curve(positive_ion_density_Combined).about.description']='Distribution of positive ions confined within the nanoparticle surfaces'
io['output.curve(positive_ion_density_Combined).xaxis.label']='z'
io['output.curve(positive_ion_density_Combined).xaxis.description']='Distance measure between the two Surfaces (z = 0 is the midpoint)'
io['output.curve(positive_ion_density_Combined).xaxis.units']='nm'
io['output.curve(positive_ion_density_Combined).yaxis.label']='Density'
io['output.curve(positive_ion_density_Combined).yaxis.description']='Density distribution of ions'
io['output.curve(positive_ion_density_Combined).yaxis.units']='M'
io['output.curve(positive_ion_density_Combined).component.xy']=(xList, yList)

io['output.curve(negative_ion_density_Combined).about.label']='Density of Negative Ions'
io['output.curve(negative_ion_density_Combined).about.group']='Density of Positive/Negative Ions'
io['output.curve(negative_ion_density_Combined).about.description']='Distribution of negative ions confined within the nanoparticle surfaces'
io['output.curve(negative_ion_density_Combined).xaxis.label']='z'
io['output.curve(negative_ion_density_Combined).xaxis.description']='Distance measure between the two Surfaces (z = 0 is the midpoint)'
io['output.curve(negative_ion_density_Combined).xaxis.units']='nm'
io['output.curve(negative_ion_density_Combined).yaxis.label']='Density'
io['output.curve(negative_ion_density_Combined).yaxis.description']='Density distribution of ions'
io['output.curve(negative_ion_density_Combined).yaxis.units']='M'
io['output.curve(negative_ion_density_Combined).component.xy']=(xListN, yListN)

		
# Close the input file handler
io.close()

sys.exit(0)

			
