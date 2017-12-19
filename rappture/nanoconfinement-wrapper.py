# ----------------------------------------------------------------------
#  Python wrapper for Nanoconfinement code.
#
#  This wrapper script calls Rappture API to fetch input values and
#   passes them to nanoconfinement tool. The wrapper can be extended
#   to do higher level data manipulation to marshall user inputs.
#
# ----------------------------------------------------------------------
import Rappture
import sys, os, commands, string

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
print " confinement_length is %s" % confinement_length

ion_diameter = io['input.group(physical).number(ion_diameter).current'].value
print " ion_diameter is %s" % ion_diameter

simulation_steps = io['input.group(computing).integer(simulation_steps).current'].value
print "simulation_steps is %s" % simulation_steps


try:
     exitStatus,stdOutput,stdError = Rappture.tools.executeCommand(
        ['md_simulation_confined_ions', '-Z', confinement_length, '-p', positive_valency, '-n', negative_valency, '-c',
         salt_concentration, '-d', ion_diameter, '-S', simulation_steps], streamOutput=True)

except:
    sys.stderr.write('Error during execution of md_simulation_confined_ions')
    sys.exit(1);
	

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

fid = open('data/p_density_profile_3.00_1_-1_0.50_0.714_5000.dat','r')
info = fid.readlines()
fid.close()

# add density profile to xy data
xList = []
yList = []
for line in info:
	proLine=" ".join(line.split())
	print proLine
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

fid = open('data/n_density_profile_3.00_1_-1_0.50_0.714_5000.dat','r')
info = fid.readlines()
fid.close()

# add density profile to xy data
xList = []
yList = []
for line in info:
	proLine=" ".join(line.split())
	print proLine
	d,m,e = proLine.split()
	xList.append(float(d))
	yList.append(float(m))

io['output.curve(negative_ion_density).component.xy']=(xList, yList)

		
# Close the input file handler
io.close()

sys.exit(0)
