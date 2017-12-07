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

try:
  out = commands.getoutput('../bin/md_simulation_confined_ions')
  print out
except:
  sys.stderr.write('Error during execution of md_simulation_confined_ions')
  exit(1);


# driver.put('output.log',out)
#
# fid = open('out.dat','r')
# info = fid.readlines()
# fid.close()
#
# # Label output graph with title, x-axis label,
# # y-axis lable, and y-axis units
# driver.put('output.curve(f12).about.label','Fermi-Dirac Factor')
# driver.put('output.curve(f12).xaxis.label','Fermi-Dirac Factor')
# driver.put('output.curve(f12).yaxis.label','Energy')
# driver.put('output.curve(f12).yaxis.units','eV')
#
# # skip over the first 4 header lines
# for line in info[6:]:
#     f,E = string.split(line[:-1])
#     f,E = float(f),float(E)
#     xy = "%g %g\n" % (f,E)
#     driver.put('output.curve(f12).component.xy',xy,append=1)
#
# os.remove('indeck'); os.remove('out.dat')
#
# Rappture.result(driver)

#Close the input file handler
user_inputs.close()
sys.exit(0)
