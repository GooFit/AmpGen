import sys
import re
file = open( sys.argv[1], "r" )
output = open( sys.argv[2], "w")
for line in file: 
    regex = "@f\$([^@]*)@f\$"
    matches = re.findall(regex, line)
    counter=0
    for match in matches: 
        replacement = "<img src=\"https://render.githubusercontent.com/render/math?math="+match +"\">"
        replacement = replacement.replace("+", "%2B") 
        line = line.replace( "@f$" + match + "@f$", replacement) 
    output.write(line)

output.close()
