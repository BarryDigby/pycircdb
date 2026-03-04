from . import config
import logging 

# Initialise the logger
logger = logging.getLogger(__name__)

def ingest(**kwargs):

    if circRNA_algorithm is not None:
        uniq = circRNA_algorithm.split(',')
        invalid = [x for x in uniq if x not in circrna_algorithms]
        if len(invalid) > 0:
            logger.critical("Invalid parameter(s) provided to `--circRNA-algorithms`: {} ".format(','.join(invalid)))
            logger.critical("Please select from the following: {}".format(','.join(circrna_algorithms)))
            return {"report": None, "config": None, "sys_exit_code": 1}

    
if __name__ == "__main__":
    ingest()
