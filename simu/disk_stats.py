import os.path
import sys

template='''
# Execute job in the queue "std.q" unless you have special requirements.
##$ -q std.q
##$ -q all_24.q

# The SGE batch system uses the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batch system assumes to find the executable in this directory.
#$ -cwd

# Redirect output stream to this file.
##$ -o output.dat

# Redirect error stream to this file.
##$ -e error.dat

# Send status information to this email address.
##$ -M oleksandr.frei@gmail.com

# Send an e-mail when the job is done.
##$ -m e

#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
##$ -pe dmp4 16

#$ -l h=mmil-compute-{}-{}.local

echo hostname
df -h /tmp
df -h /scratch

#ls /tmp

ls /scratch/simu*

rm /scratch/simu_*

'''

nodes = {5:[x for x in range(19) if (x != 8)], 7:list(range(9)), 8:list(range(1,20)), 6:list(range(9)) }

for node in nodes[int(sys.argv[1])]:
    print(sys.argv[1], node)
    with open('disk_stats.{}.{}.sh'.format(sys.argv[1], node), 'w') as f: f.write(template.format(int(sys.argv[1]), node))
    os.system('qsub disk_stats.{}.{}.sh'.format(sys.argv[1], node))
