import subprocess

def run_cmd( cmd, log ):
    print( cmd )
    subprocess.call( cmd, shell = True )
    log.write("%s\n" % ( cmd ))

