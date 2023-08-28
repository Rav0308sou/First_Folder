from _NetworkData import *
import cmath


def Polar2Rec(Vm, delta):

    V = cmath.rect(Vm, delta)
    return V

# =============================================================================
# Foramtion Of Array of Load and Npq bus:
# =============================================================================


def Npq_code(net):

    net.no_pq = net.no_bus-1-net.no_pv

    for i in range(net.no_bus):

        m = net.bus[i].code

        if net.bus[i].type != 1 and net.bus[i].type != 2:
            net.Npq.append(m)


def N1_code(net):

    for i in range(net.no_bus):

        m = net.bus[i].code

        if net.bus[i].type != 1:
            net.N1.append(m)
# =============================================================================
# Voltage Initialization:
# =============================================================================


def Voltage_Initialization(net):

    for i in range(net.no_pq):

        m = net.Npq[i]

        net.bus[m].Vm = 1.0
        net.bus[m].delta = 0.0
        V = Polar2Rec(net.bus[m].Vm, net.bus[m].delta)
        net.bus[m].V = V

    for i in range(net.no_pv):

        m = net.pv_bus[i].code

        net.bus[m].Vm = net.pv_bus[i].Vsp
        net.bus[m].delta = 0.0

        V = Polar2Rec(net.bus[m].Vm, net.bus[m].delta)
        net.bus[m].V = V

    m = net.slack_bus.code

    net.bus[m].Vm = net.slack_bus.Vsp
    net.bus[m].delta = 0.0
    V = Polar2Rec(net.bus[m].Vm, net.bus[m].delta)
    net.bus[m].V = V
# =============================================================================
# Admittance Matrix Formation:
# =============================================================================


def Admittance_Matrix(net):

    net.Y = np.zeros((net.no_bus, net.no_bus), dtype=complex)
    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        net.Y[n, n] = net.Y[n, n]+net.fd[i].y+(net.fd[i].ys/2)
        net.Y[m, m] = net.Y[m, m]+net.fd[i].y+(net.fd[i].ys/2)
        net.Y[m, n] = net.Y[m, n]-net.fd[i].y
        net.Y[n, m] = net.Y[n, m]-net.fd[i].y

    if net.no_tr > 0:

        for i in range(net.no_tr):

            m = net.tr[i].end_bus[0]
            n = net.tr[i].end_bus[1]

            a1 = net.tr[i].alpha

            net.Y[m, m] = net.Y[m, m]+((a1)**2)*(net.tr[i].y)
            net.Y[n, n] = net.Y[n, n]+net.tr[i].y
            net.Y[m, n] = net.Y[m, n]-np.conj(a1)*(net.tr[i].y)
            net.Y[n, m] = net.Y[n, m]-(a1)*(net.tr[i].y)

    if net.no_shunt > 0:

        for i in range(net.no_shunt):

            m = net.shunt[i].code
            net.Y[m, m] = net.Y[m, m]+net.shunt[i].y

    net.G = np.real(net.Y)
    net.B = np.imag(net.Y)

# =============================================================================
# Count Of Load
# =============================================================================


def Count_Load(net):

    net.no_load = 0
    for i in range(net.no_bus):
        if net.bus[i].load_status == 1:
            net.no_load += 1

# =============================================================================
# Calculation Of active Power Injection and Reactive Power injection:
# =============================================================================


def Calculation_P(net):

    for i in range(net.no_bus):
        p = 0.0
        for k in range(net.no_bus):

            p = p+net.bus[i].Vm*net.bus[k].Vm*(net.G[i, k]*np.cos(
                net.bus[i].delta-net.bus[k].delta)+net.B[i, k]*np.sin(net.bus[i].delta-net.bus[k].delta))

        net.bus[i].P = p
        # print(p)


def Calculation_Q(net):

    for i in range(net.no_bus):
        q = 0.0
        for k in range(net.no_bus):
            q = q+net.bus[i].Vm*net.bus[k].Vm*(net.G[i, k]*np.sin(
                net.bus[i].delta-net.bus[k].delta)-net.B[i, k]*np.cos(net.bus[i].delta-net.bus[k].delta))

        net.bus[i].Q = q
# =============================================================================
# Calculation Of Specified P nd Q
# =============================================================================


def Calculation_Psp(net):
    
    for i in range(net.no_bus):
        net.bus[i].Psp = net.bus[i].Pg-net.bus[i].Pd


def Calculation_Qsp(net):
    
    for i in range(net.no_bus):

        net.bus[i].Qsp = net.bus[i].Qg-net.bus[i].Qd

# =============================================================================
# deviation Of Active Power and Reactive Power Calculation:
# =============================================================================

# def Calculation_delP(net):
#     net.delP=np.zeros((net.no_bus-1,1),dtype=float)

#     for i in range(net.no_bus-1):
#         m=net.N1[i]
#         net.bus[m].delP=net.bus[m].Psp-net.bus[m].P

#         net.delP[i,0]=net.bus[m].delP


def Calculation_delP(net):

    net.delP = np.zeros((net.no_bus, 1), dtype=float)
    net.delP1 = np.zeros((net.no_bus-1, 1), dtype=float)   # reduce delP

    for i in range(net.no_bus):
        # m=net.N1[i]
        net.bus[i].delP = net.bus[i].Psp-net.bus[i].P
        net.delP[i] = net.bus[i].delP

    # for reduce delP
    for i in range(net.no_bus-1):

        m = net.N1[i]

        net.delP1[i, 0] = net.bus[m].delP


# def Calculation_delQ(net):

#     net.delQ=np.zeros((net.no_pq,1),dtype=float)

#     for i in range(net.no_pq):
#         m=net.Npq[i]
#         net.bus[m].delQ=net.bus[m].Qsp-net.bus[m].Q
#         net.delQ[i,0]=net.bus[m].delQ

def Calculation_delQ(net):

    net.delQ = np.zeros((net.no_bus, 1), dtype=float)
    net.delQ1 = np.zeros((net.no_pq, 1), dtype=float)   # reduce delQ

    for i in range(net.no_bus):
        net.bus[i].delQ = net.bus[i].Qsp-net.bus[i].Q
        net.delQ[i, 0] = net.bus[i].delQ

    for i in range(net.no_pq):
        m = net.Npq[i]
        # net.bus[m].delQ=net.bus[m].Qsp-net.bus[m].Q
        net.delQ1[i, 0] = net.bus[m].delQ

# =============================================================================
# Formation Of reduced delP and delQ:
# =============================================================================


def Foramtion_delPQ(net):

    net.delPQ = np.concatenate((net.delP1, net.delQ1), axis=0)

# =============================================================================
# Jacobian Matrix Formation:
# =============================================================================


def Jacobian_Matrix(net):

    # ========================= H-Matrix formation=================================
    net.H = np.zeros((net.no_bus-1, net.no_bus-1), dtype=float)

    for i in range(net.no_bus-1):
        n = net.N1[i]
        for k in range(net.no_bus-1):
            m = net.N1[k]
            if n == m:
                net.H[i, k] = -net.bus[n].Q-(net.bus[n].Vm**2)*net.B[n, n]
            else:
                net.H[i, k] = net.bus[n].Vm*net.bus[m].Vm*(net.G[n, m]*np.sin(
                    net.bus[n].delta-net.bus[m].delta)-net.B[n, m]*np.cos(net.bus[n].delta-net.bus[m].delta))

# ============================ N-Matrix formation =============================
    net.N = np.zeros((net.no_bus-1, net.no_pq), dtype=float)

    for i in range(net.no_bus-1):

        for k in range(net.no_pq):
            n = net.N1[i]
            m = net.Npq[k]
            if n == m:
                net.N[i, k] = net.bus[n].P+net.bus[n].Vm**2*net.G[n, n]
            else:
                net.N[i, k] = net.bus[n].Vm*net.bus[m].Vm*(net.G[n, m]*np.cos(
                    net.bus[n].delta-net.bus[m].delta)+net.B[n, m]*np.sin(net.bus[n].delta-net.bus[m].delta))

# ============================ J-Matrix formation =============================
    net.J = np.zeros((net.no_pq, net.no_bus-1), dtype=float)

    for i in range(net.no_pq):
        for k in range(net.no_bus-1):
            n = net.Npq[i]
            m = net.N1[k]
            if n == m:
                net.J[i, k] = net.bus[n].P-net.bus[n].Vm**2*net.G[n, n]

            else:
                net.J[i, k] = net.bus[n].Vm*net.bus[m].Vm*(-net.G[n, m]*np.cos(
                    net.bus[n].delta-net.bus[m].delta)-net.B[n, m]*np.sin(net.bus[n].delta-net.bus[m].delta))

# =========================== L-Matirix Formation =============================
    net.L = np.zeros((net.no_pq, net.no_pq), dtype=float)

    for i in range(net.no_pq):
        for k in range(net.no_pq):
            n = net.Npq[i]
            m = net.Npq[k]

            if n == m:
                net.L[i, k] = net.bus[n].Q-(net.bus[n].Vm**2)*net.B[n, n]
            else:
                net.L[i, k] = net.bus[n].Vm*net.bus[m].Vm*(net.G[n, m]*np.sin(
                    net.bus[n].delta-net.bus[m].delta)-net.B[n, m]*np.cos(net.bus[n].delta-net.bus[m].delta))


# ================Jacobian Matrix formation ===================================

    H_N = np.concatenate((net.H, net.N), axis=1)  # H and N matrix marged
    J_L = np.concatenate((net.J, net.L), axis=1)  # J and L matrix marged

    # Marged of H, N,J nad L matrix
    net.Jaco = np.concatenate((H_N, J_L), axis=0)


# =============================================================================
# Calculation Deviation Of angel and Voltage:
# =============================================================================

def del_DeltaVoltage(net):

    net.inv_Jaco = np.linalg.inv(net.Jaco)

    net.dDeltaV = np.dot(net.inv_Jaco, net.delPQ)

# =============================================================================
# Variables Updation :
# =============================================================================


def Updation_variables(net):

    net.dDelta = np.zeros((net.no_bus-1, 1), dtype=float)
    net.dV = np.zeros((net.no_pq, 1), dtype=float)

    for i in range(net.no_bus-1):

        net.dDelta[i, 0] = net.dDeltaV[i]

        m = net.N1[i]
        net.bus[m].delta = net.bus[m].delta+net.dDelta[i, 0]

    for i in range(net.no_pq):
        net.dV[i] = net.dDeltaV[net.no_bus-1+i]

        m = net.Npq[i]
        net.bus[m].Vm = net.bus[m].Vm+net.dV[i]

    for i in range(net.no_bus):
        V = Polar2Rec(net.bus[i].Vm, net.bus[i].delta)
        net.bus[i].V = V
# =============================================================================
# Active and Reactive Power Calculation for slack bus:
# =============================================================================


def calculation_PQ(net):

    m = net.slack_bus.code
    p = 0
    for k in range(net.no_bus):
        p = p+net.bus[m].Vm*net.bus[k].Vm*(net.G[m, k]*np.cos(
            net.bus[m].delta-net.bus[k].delta)+net.B[m, k]*np.sin(net.bus[m].delta-net.bus[k].delta))

    net.bus[m].P = p

    q = 0
    for k in range(net.no_bus):
        q = q+net.bus[m].Vm*net.bus[k].Vm*(net.G[m, k]*np.sin(
            net.bus[m].delta-net.bus[k].delta)-net.B[m, k]*np.cos(net.bus[m].delta-net.bus[k].delta))

    net.bus[m].Q = q

# =============================================================================
# Newton-Raphson Method:
# =============================================================================


def NRLF_method(net):

    net.e = 1e-6
    net.rrr = -1

    for itr in range(10):

        net.rrr = net.rrr+1

        Calculation_P(net)

        Calculation_Q(net)

        Calculation_delP(net)

        Calculation_delQ(net)

        Foramtion_delPQ(net)

        if np.max(np.abs(net.delPQ)) > net.e:

            Jacobian_Matrix(net)

            del_DeltaVoltage(net)

            Updation_variables(net)

        else:

            calculation_PQ(net)

            print(f'\n NRLF Convereged in {itr} iterations\n')

            break

# =============================================================================
# Definations of FDLF started:
# =============================================================================


def Yd_matrix(net):

    net.Yd = np.zeros((net.no_bus, net.no_bus), dtype=complex)

    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        net.Yd[n, n] = net.Yd[n, n]+1/np.complex(0, net.fd[i].x)
        net.Yd[m, m] = net.Yd[m, m]+1/np.complex(0, net.fd[i].x)
        net.Yd[n, m] = net.Yd[n, m]-1/np.complex(0, net.fd[i].x)
        net.Yd[m, n] = net.Yd[m, n]-1/np.complex(0, net.fd[i].x)

    if net.no_tr > 0:

        for i in range(net.no_tr):

            m = net.tr[i].end_bus[0]
            n = net.tr[i].end_bus[1]

            a1 = net.tr[i].alpha

            net.Yd[m, m] = net.Yd[m, m]+((a1)**2)*(net.tr[i].y)
            net.Yd[n, n] = net.Yd[n, n]+net.tr[i].y
            net.Yd[m, n] = net.Yd[m, n]-np.conj(a1)*(net.tr[i].y)
            net.Yd[n, m] = net.Yd[n, m]-(a1)*(net.tr[i].y)

    net.B1 = np.imag(net.Yd)
    net.G1 = np.real(net.Yd)


def Bd_matrix(net):

    # Bd matrix =all row and column will zero in connection with slack bus and
    # diagonal position of slack bus will contain some large number

    net.Bd = - copy.deepcopy(net.B1)
    # net.Bd = - copy.deepcopy(net.B)


    m = net.slack_bus.code
    net.Bd[m, m] = 1e8

    # for i in range(net.no_bus):
    #     m = net.slack_bus.code
    #     net.Bd[m, i] = 0.0
    #     net.Bd[i, m] = 0.0

    for i in range(net.no_bus):
        if net.Bd[i, i] == 0:
            net.Bd[i, i] = 1e8

    net.Bd_inv = np.linalg.inv(net.Bd)


def Bdd_matrix(net):

    net.Bdd = -copy.deepcopy(net.B)

    m = net.slack_bus.code
    net.Bdd[m, m] = 1e8

    for i in range(net.no_pv):
        m = net.pv_bus[i].code
        net.Bdd[m, m] = 1e8

    net.Bdd_inv = np.linalg.inv(net.Bdd)


def Calculate_delDelta(net):

    net.delPV = np.zeros((net.no_bus, 1), dtype=float)

    for i in range(net.no_bus):

        net.delPV[i, 0] = (net.delP[i, 0]/net.bus[i].Vm)

    net.delDelta = np.dot(net.Bd_inv, net.delPV)


def Calculate_delV(net):

    net.delQV = np.zeros((net.no_bus, 1), dtype=float)

    for i in range(net.no_bus):
        net.delQV[i, 0] = (net.delQ[i, 0]/net.bus[i].Vm)

    net.delV = np.dot(net.Bdd_inv, net.delQV)


def Updation_Angle(net):

    for i in range(net.no_bus):
        
        net.bus[i].delta = net.bus[i].delta+net.delDelta[i, 0]


def Updation_Voltage(net):

    for i in range(net.no_bus):
        net.bus[i].Vm = net.bus[i].Vm+net.delV[i, 0]

    for i in range(net.no_bus):
        V = Polar2Rec(net.bus[i].Vm, net.bus[i].delta)
        net.bus[i].V = V

    V = np.zeros((net.no_bus, 1), dtype=complex)
    for i in range(net.no_bus):

        V[i, 0] = net.bus[i].V

    # print(V)

def FDLF_Bd_Bdd_matrix(net):
    
    Yd_matrix(net)
    Bd_matrix(net)
    Bdd_matrix(net)

def FDLF_method(net):

    net.e = 1e-6
    net.rrr = -1

    for itr in range(1000):

        net.rrr = net.rrr+1

        Calculation_P(net)

        Calculation_delP(net)
        
        if np.max(np.abs(net.delP1)) > net.e:

            Calculate_delDelta(net)

            Updation_Angle(net)

            # Updation_Voltage(net)
        Calculation_Q(net)

        Calculation_delQ(net)
         
       
        if np.max(np.abs(net.delP1)) and np.max(np.abs(net.delQ1)) > net.e:
            
            Calculate_delV(net)

            Updation_Voltage(net)

        else:
            
            calculation_PQ(net)
           
            print(f'\n FDLF Convereged in {itr} iterations\n')

            break


def FDLF1_method(net):

    net.e = 1e-6
    
    Calculation_P(net)

    Calculation_delP(net)
  
    net.rrr = -1

    for itr in range(100):

        net.rrr = net.rrr+1

        if np.max(np.abs(net.delP)) > net.e:

            Calculate_delDelta(net)

            Updation_Angle(net)

            Calculation_P(net)

            Calculation_delP(net)

            # Updation_Voltage(net)

        if np.max(np.abs(net.delP)) and np.max(np.abs(net.delQ1))>net.e:

        # if np.max(np.abs(net.delQ1)) > net.e:
            print('sourav1')
            
            Calculation_Q(net)

            Calculation_delQ(net)

            Calculate_delV(net)

            Updation_Voltage(net)

            Calculation_Q(net)

            Calculation_delQ(net)

        else:

            calculation_PQ(net)
            print('sourav')
            print(f'\nConvereged in {itr} iterations\n')

            break


# =============================================================================
# Calculation of Generator Active and Reactive Power
# =============================================================================

def Calculate_Pg(net):
    for i in range(net.no_gen):
        a=net.Ngen[i]
        net.bus[a].Pg = net.bus[a].P+net.bus[a].Pd


def Calculate_Qg(net):
    for i in range(net.no_gen):
        a=net.Ngen[i]
        net.bus[a].Qg = net.bus[a].Q+net.bus[a].Qd

# =============================================================================
# Calculation of feeder current and power:
# =============================================================================


def Calculate_fd_I(net):

    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        net.fd[i].I[0] = (net.bus[n].V-net.bus[m].V) * \
            net.fd[i].y+net.bus[n].V*net.fd[i].ys/2
        
        net.fd[i].I[1] = (net.bus[m].V-net.bus[n].V) * \
            net.fd[i].y+net.bus[m].V*net.fd[i].ys/2

        net.fd[i].S[0] = net.bus[n].V*np.conj(net.fd[i].I[0])
        net.fd[i].S[1] = net.bus[m].V*np.conj(net.fd[i].I[1])


def Calculate_fd_Power(net):

    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        net.fd[i].S[0] = net.bus[n].V*np.conj(net.fd[i].I[0])
        net.fd[i].S[1] = net.bus[m].V*np.conj(net.fd[i].I[1])

        net.fd[i].P[0] = np.real(net.fd[i].S[0])
        net.fd[i].P[1] = np.real(net.fd[i].S[1])

        net.fd[i].Q[0] = np.imag(net.fd[i].S[0])
        net.fd[i].Q[1] = np.imag(net.fd[i].S[1])
        
        # for contingencey:
            
        # net.fd[i].S_C[0]=net.fd[i].S[0]
        # net.fd[i].S_C[1]=net.fd[i].S[1]
        
        # net.fd[i].P_C[0]=net.fd[i].P[0]
        # net.fd[i].P_C[0]=net.fd[i].P[0]
        
        # net.fd[i].Q_C[0]=net.fd[i].Q[0]
        # net.fd[i].Q_C[0]=net.fd[i].Q[0]

# =============================================================================
# calculation of Transformer currenr and power:
# =============================================================================


def tr_Current(net):

    for i in range(net.no_tr):

        a1 = net.tr[i].alpha

        m = net.tr[i].end_bus[0]
        n = net.tr[i].end_bus[1]

        net.tr[i].I[0] = ((net.bus[m].V)*(a1**2*(net.tr[i].y)) +
                          (net.bus[n].V)*(-np.conj(a1)*net.tr[i].y))
        net.tr[i].I[1] = ((net.bus[n].V)*(net.tr[i].y) +
                          (net.bus[m].V)*(-a1*net.tr[i].y))

        net.tr[i].S[0] = (net.bus[m].V)*np.conj(net.tr[i].I[0])
        net.tr[i].S[1] = (net.bus[n].V)*(np.conj(net.tr[i].I[1]))

        net.tr[i].P[0] = np.real(net.tr[i].S[0])
        net.tr[i].P[1] = np.real(net.tr[i].S[1])

        net.tr[i].Q[0] = np.imag(net.tr[i].S[0])
        net.tr[1].Q[1] = np.imag(net.tr[i].S[1])

# =============================================================================
# Load Code
# =============================================================================


def load_code(net):

    # net.load_S=np.zeros((net.no_load,1),dtype=complex)

    for i in range(net.no_bus):
        #
        if net.bus[i].Pd > 0 or net.bus[i].Qd > 0:

            net.Nload.append(i)
# =============================================================================
# Load Current and Power Calculation:
# =============================================================================


def load_S_I(net):

    net.load_S = np.zeros((net.no_load, 1), dtype=complex)
    net.load_I = np.zeros((net.no_load, 1), dtype=complex)

    for i in range(net.no_bus):

        if net.bus[i].Pd > 0 and net.bus[i].Qd > 0:

            net.bus[i].load_S = complex(net.bus[i].Pd, net.bus[i].Qd)
            net.bus[i].load_I = np.conj(net.bus[i].load_S/net.bus[i].V)

    for i in range(net.no_load):
        m = net.Nload[i]

        net.load_S[i, 0] = net.bus[m].load_S
        net.load_I[i, 0] = np.conj(net.load_S[i, 0]/net.bus[m].V)

# =============================================================================
# Load Impedance Calculation:
# =============================================================================


def load_y(net):

    net.load_y = np.zeros((net.no_load, 1), dtype=complex)

    for i in range(net.no_load):
        m = net.Nload[i]

        net.load_y[i, 0] = net.load_I[i, 0]/net.bus[m].V

# =============================================================================
# Injection_current Calculation:
# =============================================================================


def Injection_Current(net):

    net.Ii = np.zeros((net.no_bus, 1), dtype=complex)

    for i in range(net.no_bus):

        for k in range(net.no_bus):
            net.Ii[i, 0] = net.Ii[i, 0]+net.Y[i, k]*net.bus[k].V

        net.bus[i].Ii = net.Ii[i, 0]

# =============================================================================
# Generator Current ,Power ,Voltage, EMF, and angle : Calculation and Printing:
# =============================================================================


def Ngen(net):
    # net.Ngen=np.zeros((net.no_gen,1),dtype=float)
    for i in range(net.no_bus):
        m = net.bus[i].code
        if net.bus[i].type != 3 and net.bus[i].Pg>0:
            net.Ngen.append(m)


def gen_I(net):

    gen_I = np.zeros((net.no_bus, 1), dtype=complex)

    for i in range(net.no_bus):
        gen_I[i, 0] = net.bus[i].Ii+net.bus[i].load_I

    for i in range(net.no_gen):
        m = net.Ngen[i]
        net.gen[i].I = net.bus[m].Ii+net.bus[m].load_I
    # print(gen_I)


def gen_E(net):

    for i in range(net.no_gen):
        m = net.Ngen[i]
        net.gen[i].E = net.bus[m].V+net.gen[i].I*(1/net.gen[i].y)
        net.gen[i].Em = np.abs(net.gen[i].E)
        net.gen[i].delta = np.angle(net.gen[i].E)


def gen_V(net):
    
    for i in range(net.no_gen):
        
        m = net.Ngen[i]
        net.gen[i].V = net.bus[m].V
        net.gen[i].Vm = net.bus[m].Vm


def gen_PQ(net):

    for i in range(net.no_gen):

        net.gen[i].S = net.gen[i].V*np.conj(net.gen[i].I)
        net.gen[i].P = np.real(net.gen[i].S)
        net.gen[i].Q = np.imag(net.gen[i].S)

# =============================================================================
# Calculation Of Shunt Current :
# =============================================================================


def Shunt_Current(net):

    for i in range(net.no_shunt):

        m = net.shunt[i].code

        net.shunt[i].I = (net.bus[m].V)*(net.shunt[i].y)
        net.shunt[i].S = net.bus[m].V*(np.conj(net.shunt[i].I))

# =============================================================================
# KCL verification:
# =============================================================================


def kcl_mismatch(net):

    net.kcl = np.zeros(net.no_bus, dtype=complex)
    net.Itotal = np.zeros(net.no_bus, dtype=complex)
    z = -1
    net.error = np.zeros(net.no_bus, dtype=complex)
    # net.errorIndex=np.zeros(10,dtype=float)

    for i in range(net.no_bus):

        m = net.bus[i].code
        s = 0

        for k in range(net.no_fd):

            if m == net.fd[k].end_bus[0] or m == net.fd[k].end_bus[1]:

                if m == net.fd[k].end_bus[0]:
                    s = s+net.fd[k].I[0]
                else:
                    s = s+net.fd[k].I[1]
                    
        print(s)
        if net.no_tr > 0:

            for k in range(net.no_tr):

                if m == net.tr[k].end_bus[0] or m == net.tr[k].end_bus[1]:

                    if m == net.tr[k].end_bus[0]:
                        s = s+net.tr[k].I[0]

                    else:
                        s = s+net.tr[k].I[1]

        if net.no_shunt > 0:
            for k in range(net.no_shunt):
                if m == net.shunt[k].code:
                    s = s+net.shunt[k].I

        net.Itotal[i] = s
        print(s)
        net.kcl[i] = net.bus[m].Ii-net.Itotal[i]

        if net.kcl[i] > 0.1:
            z = z+1
            net.error[z] = net.kcl[i]
            net.errorIndex.append(i+1)

    # net.KCL_list[net.r,0]=np.max(np.abs(net.kcl))

    kcla = np.max(np.abs(net.kcl))
    # print(kcla)

# =============================================================================
# Loss Mismatch Calculation:
# =============================================================================


def Loss2_con(net):

    net.ss = 0.0

    fd_loss = np.zeros((net.no_fd, 1), dtype=complex)
    tr_loss = np.zeros((net.no_tr, 1), dtype=complex)

    for i in range(net.no_bus):
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)

    for i in range(net.no_fd):
        fd_loss[i, 0] = net.fd[i].S_C[0]+net.fd[i].S_C[1]
        net.ss += fd_loss[i]
    for i in range(net.no_tr):
        tr_loss[i, 0] = net.tr[i].S[0]+net.tr[i].S[1]
        net.ss += tr_loss[i]

    for i in range(net.no_shunt):
        net.ss += net.shunt[i].S
    ss1 = 0.0
    for i in range(net.no_bus):
        ss1 += net.bus[i].S

    print('\n')
    print('PowerLoss_Mismatch')
    print(np.abs(ss1-net.ss))

def Loss2(net):

    net.ss = 0.0

    fd_loss = np.zeros((net.no_fd, 1), dtype=complex)
    tr_loss = np.zeros((net.no_tr, 1), dtype=complex)

    for i in range(net.no_bus):
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)

    for i in range(net.no_fd):
        fd_loss[i, 0] = net.fd[i].S[0]+net.fd[i].S[1]
        net.ss += fd_loss[i]
    for i in range(net.no_tr):
        tr_loss[i, 0] = net.tr[i].S[0]+net.tr[i].S[1]
        net.ss += tr_loss[i]

    for i in range(net.no_shunt):
        net.ss += net.shunt[i].S
    ss1 = 0.0
    for i in range(net.no_bus):
        ss1 += net.bus[i].S

    print('\n')
    print('PowerLoss_Mismatch')
    print(np.abs(ss1-net.ss))


def Loss_from_InjectionPower(net):
    
    Pgen=0
    for i in range(net.no_gen):
        Pgen+=net.gen[i].P
        
    Pd=0
    
    for i in range(net.no_bus):
        Pd+=net.bus[i].Pd
    
    Pi=0.0    
    for i in range(net.no_bus):
        Pi+=net.bus[i].P
    
    Qgen=0
    for i in range(net.no_gen):
        Pgen+=net.gen[i].Q
        
    Qd=0
    
    for i in range(net.no_bus):
        Pd+=net.bus[i].Qd
    
    Qi=0.0    
    for i in range(net.no_bus):
        Qi+=net.bus[i].Q
    
    net.Si=complex(Pi,Qi)
    # print(net.Si)
    
    net.Bus_Si=complex(0,0)
    
    for i in range(net.no_bus):
        net.Bus_Si+=net.bus[i].S
    
    
    net.Loss_mismatch=net.Si-net.Bus_Si
    
# =============================================================================
# Production Cost Calculation:
# =============================================================================


def Cost_calculation1(net):

    Pg = np.zeros((net.no_gen, 1), dtype=float)
    for i in range(net.no_gen):
        # m=net.gen[i].bus
        Pg[i, 0] = net.gen[i].P

    net.cost_fdlf = 0.0
    for i in range(net.no_gen):
        net.cost_fdlf = net.cost_fdlf + \
            net.gen[i].a*(Pg[i, 0])**2+net.gen[i].b*(Pg[i, 0])+net.gen[i].c

    print('\n')
    print(f'Total Production Cost through FDLF :  {net.cost_fdlf}')


def Cost_calculation2(net):

    Pg = np.zeros((net.no_gen, 1), dtype=float)
    for i in range(net.no_gen):
        m = net.gen[i].bus
        Pg[i, 0] = net.Pg[m, 0]
        # Pg[i,0]=net.Pg[i,0]

    net.cost_opf = 0.0
    for i in range(net.no_gen):
        net.cost_opf += net.gen[i].a * \
            (Pg[i, 0])**2+net.gen[i].b*(Pg[i, 0])+net.gen[i].c

    print('\n')
    # print(f'Total Production Cost :  {net.cost1}')


def Cost_calculation3(net):

    Pg = np.zeros((net.no_gen, 1), dtype=float)
    
    for i in range(net.no_gen):
        m=net.gen[i].bus
        Pg[i,0]=net.bus[m].Pg
        # Pg[i, 0] = net.Pg[i, 0]

    net.cost_opf = 0.0
    for i in range(net.no_gen):
        net.cost_opf = net.cost_opf + \
            net.gen[i].a*(Pg[i, 0])**2+net.gen[i].b*(Pg[i, 0])+net.gen[i].c
    
    print(Pg)
    print('\n')
    # print(f'Total Production Cost :  {net.cost_opf}')


# =============================================================================
# Counting of conderser:
# =============================================================================

def Gen_BusCode_without_condenser(net):
    
    c = -1
    net.Gen_Code=[]
    net.Gen_MVA=[]

    for i in range(net.no_bus):
        if (net.bus[i].Pg != 0 and net.bus[i].type == 2) :
        # if net.bus[i].type == 2:
            c = c+1
            net.Gen_Code.append(i+1)
            # if i==net.gen[i].bus:
            # net.Gen_MVA.append(net.gen[c].rated_mva)
        elif (net.bus[i].Pg != 0 and net.bus[i].type==1):
            net.Gen_Code.append(i+1)
            c=c+1
            net.Gen_MVA.append(net.gen[c].rated_mva)
    
        
    # net.gen_serialNo=[]
    
    # for i in range(19):
        
    #    m=net.Gen_Code[i] 
    #    net.gen_serialNo.append(k)
        
                
    net.no_condenser = net.no_gen-c


# =============================================================================
# line flow calculation from line power formulea :
# =============================================================================

def Line_Flow(net):
    
    net.fd_Pij = np.zeros((net.no_fd, 1), dtype=float)

    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        g = np.real(net.fd[i].y)
        b = np.imag(net.fd[i].y)

        net.fd_Pij[i, 0] = (net.bus[n].Vm**2)*g-net.bus[n].Vm*net.bus[m].Vm*(g*np.cos(
            net.bus[n].delta-net.bus[m].delta)+b*np.sin(net.bus[n].delta-net.bus[m].delta))

# =============================================================================
# Modified Admittance Matrix formation for contingencey:
# =============================================================================

def Y_matrix_contingency(net):
    
    # net.r=0
    r=net.r
    
    # =============================================================================
    # Modified Bd matrix formation:
    # =============================================================================
        
    Yd_matrix(net)        # call the modified matrix
    
    n=net.fd[r].end_bus[0]
    m=net.fd[r].end_bus[1]
    
    a=1/(1j*net.fd[r].x)
    
    net.Yd[n,n]=net.Yd[n,n]-a
    net.Yd[m,m]=net.Yd[m,m]-a
    net.Yd[n,m]=net.Yd[n,m]+a
    net.Yd[m,n]=net.Yd[m,n]+a
    
    
    net.B1=np.imag(net.Yd)
    net.G1=np.real(net.Yd)
    
    Bd_matrix(net)
    
    # =============================================================================
    # Modified Bdd matrix formation:
    # =============================================================================
        
    Admittance_Matrix(net)   # call the original matrix again:
    
    net.Y[m,m]=net.Y[m,m]-net.fd[r].y-((net.fd[r].ys)/2)
    net.Y[n,n]=net.Y[n,n]-net.fd[r].y-((net.fd[r].ys)/2)
    net.Y[m,n]=net.Y[m,n]+net.fd[r].y
    net.Y[n,m]=net.Y[n,m]+net.fd[r].y
    
    net.G=np.real(net.Y)
    net.B=np.imag(net.Y)
    
    Bdd_matrix(net)
    
    
def Calculate_fd_I_Contingency(net):

    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        net.fd[i].I[0] = (net.bus[n].V-net.bus[m].V) * \
            net.fd[i].y+net.bus[n].V*net.fd[i].ys/2
        net.fd[i].I[1] = (net.bus[m].V-net.bus[n].V) * \
            net.fd[i].y+net.bus[m].V*net.fd[i].ys/2
            
            
        if net.r==i:
            net.fd[i].I[0]=0.0
            net.fd[i].I[1]=0.0
            
        # net.fd[i].S[0] = net.bus[n].V*np.conj(net.fd[i].I[0])
        # net.fd[i].S[1] = net.bus[m].V*np.conj(net.fd[i].I[1])
    
def Calculate_fd_Power_Contingencey(net):

    for i in range(net.no_fd):

        n = net.fd[i].end_bus[0]
        m = net.fd[i].end_bus[1]

        net.fd[i].S_C[0] = net.bus[n].V*np.conj(net.fd[i].I[0])
        net.fd[i].S_C[1] = net.bus[m].V*np.conj(net.fd[i].I[1])

        net.fd[i].P_C[0] = np.real(net.fd[i].S_C[0])
        net.fd[i].P_C[1] = np.real(net.fd[i].S_C[1])

        net.fd[i].Q_C[0] = np.imag(net.fd[i].S_C[0])
        net.fd[i].Q_C[1] = np.imag(net.fd[i].S_C[1])
        
       
def Percentage_Load(net):
    
    net.percentage_load=np.zeros((net.no_fd,1),dtype=float)
    
    for i in range(net.no_fd):
        if net.r!=i:
            
            net.percentage_load[i,0]=np.abs((net.fd[i].P[0]-net.fd[i].P_C[0])/net.fd[i].P[0])
        else:
            net.percentage_load[i,0]=0.0
            

          
def Overload_line_determination(net):
    
    net.Line_No_overload=[]
    net.Percentage=[]
    
    
    for i in range(net.no_fd):
        
        if i!=net.r:
        
            if np.abs(net.fd[i].P[0])<np.abs(net.fd[i].P_C[0]):
                net.Line_No_overload.append(i)
                net.Percentage.append(net.percentage_load[i,0]*100)
            
    a=(len(net.Percentage))
    net.per_dim.append(a)
    
    net.Line_No_overload1=np.zeros((a,1),dtype=float)
    net.Percentage1=np.zeros((a,1),dtype=float) 
    
    for i in range(a): 
        
        net.Line_No_overload1[i,0]=net.Line_No_overload[i]
        net.Percentage1[i,0]=net.Percentage[i]
        
    net.Per_Index=np.concatenate((net.Line_No_overload1,net.Percentage1),axis=1)
        

def Indexing_All(net):
    net.Index_Percent_All.append(net.Per_Index)    
    
    
    
def Bd_fv(net):
    
    net.Bd=- copy.deepcopy(net.B1)
    
    for i in range(net.no_bus):
        a=0
        for k in range(net.no_gen):
            if i==net.gen[k].bus:
                a+=(net.gen[k].P*1.5/(0.05*net.bus[i].Vm))
                
        a+=0.1*1.5*net.bus[i].Pd
        m=net.slack_bus.code
        net.Bd[i,m]=a
    
    net.Bd_inv = np.linalg.inv(net.Bd)

def Bdd_fv(net):
    
    net.Bdd = -copy.deepcopy(net.B)
    for i in range(net.no_bus):
        net.Bdd[i,i]+=2*0.05*net.bus[i].Qd
    
    net.Bdd_inv = np.linalg.inv(net.Bdd)
   
         
def del_P_fv(net):
    
    net.P_fv=np.zeros((net.no_bus,1),dtype=float)
    
    for i in range(net.no_bus):
        v=net.bus[i].Vm
        a=net.bus[i].Pd*(1+1.5*net.df)*(0.85+.1*v+.05*v**2)
        net.P_fv[i,0]=net.bus[i].P-(net.bus[i].Pg-(net.df*1.5*net.bus[i].Pg)/.05)+a
        # net.P_fv[i,0]=net.P_fv[i,0]/v
                   
    
def del_P_fv_d(net):
    
    net.P_fv_d=np.zeros((net.no_bus,1),dtype=float)
    net.delP1=np.zeros((net.no_bus-1,1),dtype=float)
    for i in range(net.no_bus):       
            v=net.bus[i].Vm
            net.P_fv_d[i,0]=net.P_fv[i,0]-net.bus[i].Pd*1.5*(0.85+.05*v**2)*net.ddf
            net.P_fv_d[i,0]/=v
    
    c=-1
    for i in range(net.no_bus):
        if net.slack_bus.code!=i:
            m=net.bus[i].code
            c+=1
            net.delP1[c,0]=net.P_fv_d[m,0]
    
    
    
def Calculate_delDelta_fv(net):

    net.delPV = np.zeros((net.no_bus, 1), dtype=float)

    
    net.delDelta = np.dot(net.Bd_inv, net.P_fv_d)

    
def Bd_fv1(net):
    
    net.Bd=- copy.deepcopy(net.B1)
    
    for i in range(net.no_bus):
        #a=0
        
        a=(net.bus[i].Pg*1.5/(0.05*net.bus[i].Vm))+0.1*1.5*net.bus[i].Pd
                
        
        m=net.slack_bus.code
        net.Bd[i,m]=a
    
    net.Bd_inv = np.linalg.inv(net.Bd)
   
def LODF(net):
    
    net.dlk=np.zeros((net.no_fd,net.no_fd),dtype=float)
    
    X=net.Bd_inv
    
    for l in range(net.no_fd):
        
        i = net.fd[l].end_bus[0]
        j = net.fd[l].end_bus[1]
        
        for k in range(net.no_fd):
             
             n = net.fd[k].end_bus[0]
             m = net.fd[k].end_bus[1]
             
             a=(X[i,n]-X[i,m]-X[j,n]+X[j,m])
             b=net.fd[k].x-(X[n,n]+X[m,m]-2*X[n,m])
             
             net.dlk[l,k]=(net.fd[k].x/net.fd[l].x)*(a/b)
              
             if l==k:
                 net.dlk[l,k]=0
    
    
    
    