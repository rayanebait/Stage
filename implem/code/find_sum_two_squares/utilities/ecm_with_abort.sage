from multiprocessing import Process, Pipe 
import time


def LBound(x,a=1/2,c=1):
    return exp(c*( (ln(x)**a) * (ln(ln(x))**(1-a)) ))

"""
def is_L_smooth(factors, L):
    is_smooth=True
    for i in factors:
	if i>10*L:
	    is_smooth=False
	    break
    return is_smooth
"""

def fac(recv_cand, send_fac):
    while True:
	    c = recv_cand.recv()
	    c = abs(c)

	    factors = ecm.factor(c,None)
	    send_fac.send(factors)

def ecm_with_abort(n, t = 0.5):
	recv_candidate,send_candidate = Pipe(duplex = False)
	recv_fac, send_fac = Pipe(duplex = False)
	proc = Process(target=fac, args=(recv_candidate, send_fac))
	proc.start()

	send_candidate.send(n)

	time.sleep(float(t))
	if recv_fac.poll():
	    factors = recv_fac.recv()
	else:
	    factors = None
	proc.terminate()
	return factors


if __name__ == "__main__":
    print(ecm_with_abort(2**10*3**2))
