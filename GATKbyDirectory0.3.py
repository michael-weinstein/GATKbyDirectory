#!/usr/bin/env python3
import os
import re

def checkargs():  #subroutine for validating commandline arguments
    import argparse #loads the required library for reading the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument ("-T", "--analysis_type", help = "Which tool to run?")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("-e", "--standard_min_confidence_threshold_for_emitting", help = "Minimum confidence for emitting.")
    parser.add_argument ("-c", "--standard_min_confidence_threshold_for_calling", help = "Minimum confidence for calling?")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("-R", "--reference_sequence", help = "Reference Genome?")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("-D", "--directory", help = "BAM File Directory?")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("--dontUseSoftClippedBases", help = "An option that is useful for RNAseq based variant calling", action = 'store_true')
    parser.add_argument ("-dcov", "--downsample_to_coverage", help = "Target coverage threshold for downsampling to coverage.")
    parser.add_argument ("-sct", "--scatter_count", help = "scatter")
    parser.add_argument ("-9", "--clobber", help = "Overwrite all the things without asking first.", action = 'store_true')
    parser.add_argument ("-L", "--intervals", help = "One or more genomic intervals over which to operate")
    parser.add_argument ("--gvcf", help = "gvcf mode")
    parser.add_argument ("-ls", "--liquid_schwartz", help = "Takes lots of nodes.  Runs jobs a ludicrous speeds.", action = 'store_true')
    parser.add_argument ("--email", help = "email address for sending notifications in liquid schwartz mode")
    args = parser.parse_args()  #puts the arguments into the args object
    
    if not os.path.isfile(args.reference_sequence):
        usage('Reference sequence file specified was not found.')
    if not os.path.isdir(args.directory):
        usage('Target directory specified was not found.')
    if args.intervals and not os.path.isfile(args.intervals):
        usage('Intervals file specified was not found')
    if not args.standard_min_confidence_threshold_for_emitting:  #if not given use the default value
        args.standard_min_confidence_threshold_for_emitting = "10"
    else:
        try:
            int(args.standard_min_confidence_threshold_for_emitting)
        except:
            usage('standard_min_confidence_threshold_for_emitting value given was not integer type')
    if not args.standard_min_confidence_threshold_for_calling: #if not given, use the default value
        args.standard_min_confidence_threshold_for_calling = "30"
    else:
        try:
            int(args.standard_min_confidence_threshold_for_calling)
        except:
            usage('standard_min_confidence_threshold_for_calling value given was not integer type')
    if not args.scatter_count:
        args.scatter_count = "20"
    else:
        try:
            int(args.scatter_count)
        except:
            usage('scatter_count must be an integer')
    if args.email:
        args.email = "-M " + args.email + " -m bea "
    else:
        args.email = ""
    return (args)

def usage(sin):  #This subroutine prints directions
    print ('Error: ' + sin)
    print ('This program is designed to automate the generation and submission of GATK scatter/gather jobs on the Hoffman2 cluster.')
    quit("Please correct errors and try again.")

class DirectoryItem(object):
    def __init__(self, filename,directory):
        self.filename = filename
        self.isDirectory = os.path.isdir(directory + "/" + filename)

def filelist(directory):
    import os
    allfiles = os.listdir(directory)
    items = []
    for file in allfiles:
        items.append(DirectoryItem(file,directory))
    return items

def yesanswer(question):  #asks the question passed in and returns True if the answer is yes, False if the answer is no, and keeps the user in a loop until one of those is given.  Also useful for walking students through basic logical python functions
    answer = False  #initializes the answer variable to false.  Not absolutely necessary, since it should be undefined at this point and test to false, but explicit is always better than implicit
    while not answer:  #enters the loop and stays in it until answer is equal to True
        print (question + ' (Y/N)')  #Asks the question contained in the argument passed into this subroutine
        answer = input('>>') #sets answer equal to some value input by the user
        if str(answer) == 'y' or str(answer) == 'Y':  #checks if the answer is a valid yes answer
            return True  #sends back a value of True because of the yes answer
        elif str(answer) == 'n' or str(answer) == 'N': #checks to see if the answer is a valid form of no
            return False  #sends back a value of False because it was not a yes answer
        else: #if the answer is not a value indicating a yes or no
            print ('Invalid response.')
            answer = False #set ansewr to false so the loop will continue until a satisfactory answer is given

class Job(object):
    def __init__(self, item, args):
        self.shortname = item.filename
        self.fullname = (args.directory + "/" + item.filename)
        extension = self.shortname[len(self.shortname)-4:len(self.shortname)]
        if extension != ".bam":
            self.isBamFile = False
        else:
            self.isBamFile = True
        if item.isDirectory:
            subdirectorylist = os.listdir(self.fullname)
            subdirectorybams = []
            for subdirectoryitem in subdirectorylist:
                if re.match('.*\.bam$', subdirectoryitem):
                    self.isBamFile = True
                    newFile = (self.fullname + "/" + subdirectoryitem).replace('//','/')
                    subdirectorybams.append(newFile)
            if self.isBamFile:
                writelist = True
                if os.path.isfile(self.shortname + ".list") and not args.clobber:
                    writelist = yesanswer("List file for " + item.filename + " directory already exists.  Overwrite?")
                if writelist or args.clobber:
                    listfile = open(self.shortname + ".list",'w')
                    for subdirectoryitem in subdirectorybams:
                        baifile = subdirectoryitem + ".bai"
                        altbaifile = subdirectoryitem[:-4] + ".bai"
                        if not (os.path.isfile(baifile) or os.path.isfile(altbaifile)):
                            if not yesanswer("Neither " + baifile + " nor " + altbaifile + " (BAM Index Files) were found.  Use anyway?"):
                                if not yesanswer("Do you want to continue this run?  No jobs have been sent to qsub yet."):
                                    quit("Goodbye!")
                                else:
                                    continue
                        listfile.write(os.getcwd() + "/" + subdirectoryitem + "\n")
                    listfile.close()
                    self.fullname = self.shortname + ".list"
                else:
                    self.isBamFile = yesanswer("Include existing list file for " + item.filename + " in this run?")
        if self.isBamFile:
            if os.path.isfile(self.fullname + ".vcf") and not args.clobber:
                self.isBamFile = yesanswer(self.shortname + ".vcf already exists.  Overwrite?")
                if not self.isBamFile:
                    if not yesanswer("Do you want to continue this run?  No jobs have been sent to qsub yet."):
                        quit("Goodbye!")
            if not os.path.isfile(self.fullname + ".bai") and not os.path.isfile(self.fullname[0:-4] + ".bai") and not self.fullname[len(self.fullname)-5:len(self.fullname)] == ".list":
                self.isBamFile = yesanswer(self.shortname + ".bai (BAM Index File) is missing.  Submit job anyway?")
                if not self.isBamFile:
                    if not yesanswer("Do you want to continue this run?  No jobs have been sent to qsub yet."):
                        quit("Goodbye!")
      
    def create (self, args): #create a the scatter/gather scala object
        # self.cleanfilename = re.sub('\W','',self.shortname)
        # self.cleanfilename = re.sub('_','',self.cleanfilename)
        if os.path.isfile(self.shortname + ".scatter.scala") and not args.clobber:
            if not yesanswer(self.shortname + ".scatter.scala already exists.  Overwrite?"):
                return False
        scala = open(self.shortname + ".scatter.scala",'w')
        scala.write("import org.broadinstitute.gatk.queue.QScript" + "\n")
        scala.write("import org.broadinstitute.gatk.queue.extensions.gatk._" + "\n")
        
        scala.write("class callVariants extends QScript {" + "\n")
        scala.write("\tdef script() {"  + "\n")
        if (args.analysis_type).lower() == 'haplotypecaller':
            analysistype = "HaplotypeCaller"
            objectname = "hc"
        elif (args.analysis_type).lower() == 'unifiedgenotyper':
            analysistype = "UnifiedGenotyper"
            objectname = "genotyper"
        else:
            quit("Invalid Analysis type.  Must be either UnifiedGenotyper or HaplotypeCaller." + "\n")
        scala.write("\t\tval " + objectname + " = new " + analysistype + "\n")
        scala.write("\t\t" + objectname + ".reference_sequence = new File (\"" + args.reference_sequence + "\")" + "\n")
        if args.standard_min_confidence_threshold_for_emitting != "":
            scala.write("\t\t" + objectname + ".standard_min_confidence_threshold_for_emitting = " + args.standard_min_confidence_threshold_for_emitting + "\n")
        if args.standard_min_confidence_threshold_for_calling != "":
            scala.write("\t\t" + objectname + ".standard_min_confidence_threshold_for_calling = " + args.standard_min_confidence_threshold_for_calling + "\n")
        scala.write("\t\t" + objectname + ".input_file :+= new File (\"" + self.fullname + "\")" + "\n")
        scala.write("\t\t" + objectname + ".out = new File (\"" + self.fullname + ".vcf\")" + "\n")
        if args.dontUseSoftClippedBases:
            scala.write("\t\t" + objectname + ".dontUseSoftClippedBases = true" + "\n")
        # if args.number_cpu_threads_per_data_thread:
        #     scala.write("\t\t" + objectname + ".number_cpu_threads_per_data_thread = " + args.number_cpu_threads_per_data_thread + "\n")
        if args.downsample_to_coverage:
            scala.write("\t\t" + objectname + ".downsample_to_coverage = " + args.downsample_to_coverage + "\n")
        if args.intervals:
            scala.write("\t\t" + objectname + ".intervals :+= new File (\"" + args.intervals + "\")" + "\n")
        scala.write("\t\t" + objectname + ".scatterCount = "+ args.scatter_count + "\n")
        scala.write("\t\t" + objectname + ".memoryLimit = 2" + "\n")
        scala.write("\t\t" + "add(" + objectname + ")" + "\n")
        scala.write("\t" + "}" + "\n")
        scala.write("}" + "\n")
        scala.close()
        return True
    
    def execute (self, args): #sends the scatter/gather job to qsub
        command = 'java -Xmx1g -Djava.io.tmpdir=tmp -jar /u/local/apps/gatk-queue/3.2.2/Queue.jar -S ' + self.shortname +'.scatter.scala -startFromScratch -qsub -jobResReq "h_data=4g,h_rt=24:00:00" -run'
        if not args.liquid_schwartz:
            os.system(command)
            return True
        else:
            bashFileName = self.shortname + 'GATKBDsubmission.sh'
            bashFile = open(bashFileName, 'w')
            bashFile.write('/bin/bash\n')
            bashFile.write(command)
            bashFile.close()
            qsubmission = "qsub -cwd -V -N " + self.shortname + " -l h_data=4G,time=24:00:00 -pe shared 2 " + args.email + bashFileName
            os.system(qsubmission)

def main():
    print ("Checking command line arguments...", end = "")
    args = checkargs()
    # if args.clobber:
    #     if not yesanswer('Clobber set, files may be overwritten without asking.  Continue?'):
    #         quit("Goodbye!")
    print ("OK\nGetting target directory contents...", end = "")
    directorycontents = filelist(args.directory) #returns an array of objects with filename and if it is a directory
    print ("OK\nCreating a list of GATK jobs...", end = "")
    jobs = []
    for item in directorycontents:
        jobs.append(Job(item, args))
    print ("OK\nCreating scatter/gather scala objects...")
    for job in jobs:
        if job.isBamFile:
            print("for " + job.shortname)
            job.create(args)
    print ("OK\nSubmitting jobs to queue...")
    for job in jobs:
        if job.isBamFile:
            print("\n\n\nfor " + job.shortname)
            job.execute(args)
    print ("OK")
    quit("Done!")

main()