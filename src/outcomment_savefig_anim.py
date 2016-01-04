import os
path = os.getcwd()
print 'path: ', path
for chap in os.listdir(path):
    print "dir/file: ", chap
    if os.path.isdir(chap):
        for file_name in os.listdir(chap):
            if file_name[-3:]== '.py':
                old_file = path+'/'+chap+'/'+file_name
                new_file = path+'/'+chap+'/'+file_name + '.tmp'
                f1 = open(old_file, 'r')
                f2 = open(new_file, 'w')
                linenumber = 0
                for line in f1:
                    linenumber +=1
                    if ('savefig(' in line) and (line[0] != '#'):
                        print "edditing file: '{0} on linenumber {1}. Out-Comment line".format(old_file, linenumber)
                        f2.write('#'+line)
                    elif ('Writer = animation' in line) and (line[0] != '#'):
                        print "edditing file: '{0} on linenumber {1}. Out-Comment line".format(old_file, linenumber)
                        f2.write('#'+line)
                    elif ('writer = Writer' in line) and (line[0] != '#'):
                        print "edditing file: '{0} on linenumber {1}. Out-Comment line".format(old_file, linenumber)
                        f2.write('#'+line)
                    elif ('anim.save' in line) and (line[0] != '#'):
                        print "edditing file: '{0} on linenumber {1}. Out-Comment line".format(old_file, linenumber)
                        f2.write('#'+line)
                    else:
                        f2.write(line)
                f1.close()
                f2.close()
                os.remove(old_file)
                os.rename(new_file, old_file)
