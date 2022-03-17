#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao

import time
import Sniper

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Run Tao Detector Simulation Analysis.')
    parser.add_argument("--evtmax", type=int, default=10, help='events to be processed')
    parser.add_argument("--seed", type=int, default=42, help='seed')
    parser.add_argument("--temp_radius", type = int, default = 0, help= "Charge Template Radius")
    parser.add_argument("--charge_template_file", default = "Alpha_charge_template", help='close dark noise')
    parser.add_argument("--close_dark_noise", action = "store_true", help='close dark noise')
    parser.add_argument("--close_inter_ct", action = "store_true", help='close internal cross talk')
    parser.add_argument("--close_charge_resolution", action = "store_true", help='close charge resolution')
    parser.add_argument("--input", default=None, nargs='+', help="specify input filename")
    parser.add_argument("--output", default="user-ana-output.root", help="specify output filename")
    parser.add_argument("--algorithm", default="ChargeTemplateRec", help="specify the analysis algorithm")

    return parser

if __name__ == "__main__":

    t0= time.clock()

    parser = get_parser()
    args = parser.parse_args()
    print args

    task = Sniper.Task("task")
    task.setEvtMax(args.evtmax)
    #task.setLogLevel(0)

    # = random svc =
    import RandomSvc
    rndm = task.createSvc("RandomSvc")
    rndm.property("Seed").set(args.seed)

    # = rootio =
    import RootIOSvc
    if not args.input:
        print("Please provide an input file for analysis")
        exit(-1)
    ris = task.createSvc("RootInputSvc/InputSvc")
    ris.property("InputFile").set(args.input)

    # = BufferMemMgr =
    import BufferMemMgr
    bufMgr = task.createSvc("BufferMemMgr")
    bufMgr.property("TimeWindow").set([0, 0]);

    # = geometry service =
    import Geometry
    simgeomsvc = task.createSvc("SimGeomSvc")

    # = root writer =
    import RootWriter
    rootwriter = task.createSvc("RootWriter")
    rootwriter.property("Output").set({"ANASIMEVT": args.output})

    # = ana detsim =
    # Sniper.loadDll("libTaoSimAna.so")
    Sniper.loadDll("libChargeTemplatePos.so")
    anadetsimalg = task.createAlg(args.algorithm)
    anadetsimalg.property("CloseInterCT").set(args.close_inter_ct)
    anadetsimalg.property("CloseChargeResolution").set(args.close_charge_resolution)
    if args.algorithm == "MakeChargeTemplate":
        anadetsimalg.property("TempRadius").set(args.temp_radius)
        anadetsimalg.property("CloseDarkNoise").set(True)
    else:
        anadetsimalg.property("CloseDarkNoise").set(args.close_dark_noise)
        anadetsimalg.property("ChargeTemplateFile").set(args.charge_template_file)

    task.show()
    task.run()

    t1 = time.clock() - t0
    print "==================================="
    print "run.py INFO: time elapsed: ", "{:10.4f}".format(t1), "s =", "{:10.4f}".format(t1/60), "min"
