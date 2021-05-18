#!/usr/bin/python3
import os


if __name__ == "__main__":

    region = "hippocampus"
    result_folder = "../exp"
    # os.system("sshpass -p katia123 rsync -avz \
    #           kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/images .")
    # os.system("sshpass -p katia123 rsync -avz \
    #           kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/atlas .")
    # os.system("sshpass -p katia123 rsync -avz \
    #           kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/scripts* .")
    # os.system("sshpass -p katia123 rsync -avz \
    #           kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/template .")
    # os.system("sshpass -p katia123 rsync -avz \
    #           kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/experiment .")
    print("sshpass -p katia123 rsync -avz \
              kpoloni@bip-server3:/databases/data2/IXI/mesh/r8 mesh/IXI")


    # os.system("sshpass -p katia123 rsync -avz \
    #           experiment kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/")
    # for fold in range(5):
    #     os.system("sshpass -p katia123 rsync -avz \
    #               kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/experiment/hippocampus/" +
    #               str(fold) + "/final_results_acc2 experiment/hippocampus/" +
    #               str(fold))

        # os.system("sshpass -p katia123 rsync -avz \
        #           kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/exp/hippocampus/" +
        #           str(fold) + "/desc_eval_images_* exp/hippocampus/" +
        #           str(fold))


        # os.system("sshpass -p katia123 rsync -avz \
        #       kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/exp/hippocampus/" +
        #           str(fold) + "/final_results_acc/* exp/hippocampus/" +
        #           str(fold) + "/final_results_acc/")
        
        # os.system("sshpass -p katia123 rsync -avz \
        #       kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/exp/hippocampus/" +
        #           str(fold) + "/final_results_acc_comp/* exp/hippocampus/" +
        #           str(fold) + "/final_results_acc_comp/")
        
        # os.system("sshpass -p katia123 rsync -avz \
        #       kpoloni@bip-server3:~/workspace/experiments/concomitant/bspline/exp/hippocampus/" +
        #           str(fold) + "/final_results_joinN/* exp/hippocampus/" +
        #           str(fold) + "/final_results_joinN/")

        

