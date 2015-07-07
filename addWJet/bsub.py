#!/usr/bin/env python
  
import os
import optparse
 
def main():
 
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
       
    parser.set_defaults(overwrite=False)
    parser.add_option('-r', '--envDir',    dest='envDir',   help='Full path to CMSSW release',)
    parser.add_option('-w', '--workDir',   dest='workDir',  help='Work directory dir',)
    parser.add_option('-i', '--inDir',   dest='inDir',  help='Work directory with input/ dir',)
    parser.add_option('-o', '--outDir',   dest='outDir',  help='Work directory with output/ dir',)
    parser.add_option('-l', '--lumi',      dest='lumi',     help='Integrated luminosity', default='19.468')
  
    (opt, args) = parser.parse_args()
  
    if opt.envDir is None:
        parser.error('No ShapeAnalysis/test directory defined')
    if opt.workDir is None:
        parser.error('No work directory with shape.py defined')
    if opt.inDir is None:
        parser.error('No work directory with shape.py defined')
    if opt.outDir is None:
        parser.error('No work directory with shape.py defined')
        
    envDir = opt.envDir
    workDir = opt.workDir
    inDir = opt.inDir
    outDir = opt.outDir


    print 'Env directory: '+envDir
    print 'Work directory: '+workDir
    print 'Input directory: '+inDir
    print 'Output directory: '+outDir

    os.system('mkdir -p '+outDir+'/jobs')

    allFiles = os.listdir(inDir)
    jobs = []
    for file in allFiles:
        
        title = 'addWJetsWeights_'+file
        print 'sub_'+title+'.sh'
        subfile = open(outDir+'/jobs/sub_'+title+'.sh','w')

        subfile.write('#!/bin/bash\n')
        subfile.write('#$ -N '+title+'.root'+'\n')
        subfile.write('#$ -q all.q\n')
        subfile.write('#$ -cwd\n')

        subfile.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
        #subfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
        subfile.write('export SCRAM_ARCH=slc5_amd64_gcc462\n')
        subfile.write('cd '+envDir+'\n')
        subfile.write('eval `scramv1 ru -sh`\n')

        subfile.write('cd '+workDir+'\n')

        #subfile.write('root -b -q runWJets.C\\(\\"'+file+'\\",'+opt.lumi+',\\"input\\",\\"output\\",true\\)') ;

        subfile.write('root -b -q runWJets.C\\(\\"'+file+'\\",'+opt.lumi+',\\"'+inDir+'\\",\\"'+outDir+'\\",true\\)') ;

        subfile.close()

        jobs.append(outDir+'/jobs/sub_'+title+'.sh')
                
## submit the jobs
    for job in jobs:
        #os.system('qsub '+job)
        os.system('chmod +x '+job)
        os.system('bsub -R "pool>30000" -q 8nh '+job)
        print('qsub '+job)

                
if __name__ == '__main__':
    main()
