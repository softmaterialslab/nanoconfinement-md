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
user_inputs = Rappture.PyXml(sys.argv[1])

salt_concentration = user_inputs['input.group(physical).number(salt_concentration).current'].value
print "salt concentration is %s" % salt_concentration

positive_valency = user_inputs['input.group(physical).integer(positive_valency).current'].value
print "positive valency is %s" % positive_valency

negative_valency = user_inputs['input.group(physical).integer(negative_valency).current'].value
print "negative_valency  s %s" % negative_valency

confinement_length = user_inputs['input.group(physical).number(confinement_length).current'].value
print " confinement_length is %s" % confinement_length

ion_diameter = user_inputs['input.group(physical).number(ion_diameter).current'].value
print " ion_diameter is %s" % ion_diameter

simulation_steps = user_inputs['input.group(computing).integer(simulation_steps).current'].value
print "simulation_steps is %s" % simulation_steps

# Close the input file handler
user_inputs.close()

try:
    print Rappture.tools.executeCommand(
        ['md_simulation_confined_ions', '-Z', confinement_length, '-p', positive_valency, '-n', negative_valency, '-c',
         salt_concentration, '-d', ion_diameter, '-S', simulation_steps], streamOutput=True)
    # exitStatus, stdOut, stdErr = Rappture.tools.executeCommand(
    #     ['md_simulation_confined_ions', '-Z', confinement_length, '-p', positive_valency, '-n', negative_valency, '-c',
    #      salt_concentration, '-d', ion_diameter, '-S', simulation_steps], streamOutput=True)
except:
    sys.stderr.write('Error during execution of md_simulation_confined_ions')
    sys.exit(1);
	


# Label output graph with title, x-axis label,
# y-axis lable, and y-axis units
driver.put('output.curve(positive_ion_density).about.label','Density of Positive Ions')
driver.put('output.curve(positive_ion_density).about.description','Distribution of positive ions confined within the nanoparticle surfaces')
driver.put('output.curve(positive_ion_density).xaxis.label','z')
driver.put('output.curve(positive_ion_density).xaxis.description','Distance measure between the two Surfaces (z = 0 is the midpoint)')
driver.put('output.curve(positive_ion_density).xaxis.units','nm')
driver.put('output.curve(positive_ion_density).yaxis.label','Density')
driver.put('output.curve(positive_ion_density).yaxis.description','Density distribution of ions')
driver.put('output.curve(positive_ion_density).yaxis.units','M')

fid = open('data/p_density_profile_3.00_1_-1_0.50_0.714_5000.dat','r')
info = fid.readlines()
fid.close()

# add density profile to xy data
for line in info:
	driver.put('output.curve(positive_ion_density).component.xy',line,append=1)
	
#os.remove('indeck'); os.remove('out.dat')
	
Rappture.result(driver)

print driver.xml()
ass
sys.exit(0)
