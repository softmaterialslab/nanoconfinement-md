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

sys.exit(0)
