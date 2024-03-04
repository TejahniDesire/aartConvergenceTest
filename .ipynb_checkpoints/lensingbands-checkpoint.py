from aart_func import *
from params import * 
import params


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('dx0', type=float)
parser.add_argument('dx1', type=float)
parser.add_argument("dx2", type=float)
args = parser.parse_args()

dx0 = args.dx0
dx1 = args.dx1
dx2 = args.dx2
params.dx0 = dx0
params.dx1 = dx1
params.dx2 = dx2

print("Computing the lensing bands")
lb.lb(dx0, dx1, dx2)