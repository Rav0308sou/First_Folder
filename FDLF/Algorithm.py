from LoadFlow_Definations import *
from Printing_Part import *
from SQP1_Definations import *
import matplotlib.pyplot as plt

plt.close("all")
print("\033[H\033[J")

# ====================================================
bus = 30
net = Network()
data_file = 'IEEE' + str(bus) + ".xlsx"
read_network_data(net,data_file)
# =====================================================

def change_load(net):
    for i in range(net.no_bus):
        net.bus[i].Pd=1.2*net.bus[i].Pd
        # net.bus[i].Qd=1.2*net.bus[i].Qd
        
# change_load(net)

Gen_BusCode_without_condenser(net)

Npq_code(net)

N1_code(net)

Ngen(net)

load_code(net)

Count_Load(net)

Voltage_Initialization(net)

Admittance_Matrix(net)

Calculation_Psp(net)

Calculation_Qsp(net)

# ===============================================

FDLF_Bd_Bdd_matrix(net)

FDLF_method(net)

# ================================================

Calculate_Pg(net)

Calculate_Qg(net)

Calculate_fd_I(net)

Calculate_fd_Power(net)

tr_Current(net)

load_S_I(net)

load_y(net)

Injection_Current(net)

gen_I(net)

gen_E(net)

gen_V(net)

gen_PQ(net)

Shunt_Current(net)

kcl_mismatch(net)

Print(net)

print(f'Total iteration requried to converged for FDLF :  {net.rrr}')

Loss2(net)    

Cost_calculation1(net)

LODF(net)

def max_powerlimit(net):
    net.fd_Pmax=np.zeros((net.no_fd,1),dtype=float)
    
    for i in range(net.no_fd):
        net.fd_Pmax[i,0]=1.2*net.fd[i].P[0]
        
    
max_powerlimit(net)

def store_voltage(net):
    net.Voltage1=np.zeros((net.no_bus,1),dtype=float)
    net.delta1=np.zeros((net.no_bus,1),dtype=float)
    net.fd_P1=np.zeros((net.no_fd,1),dtype=float)
    for i in range(net.no_bus):
        net.Voltage1[i,0]=net.bus[i].Vm
        net.delta1[i,0]=net.bus[i].delta
        
    for i in range(net.no_fd):
        net.fd_P1[i,0]=net.fd[i].P[0]
        

store_voltage(net)

def store_voltage2(net):
    net.Voltage2=np.zeros((net.no_bus,1),dtype=float)
    net.delta2=np.zeros((net.no_bus,1),dtype=float)
    net.fd_P2=np.zeros((net.no_fd,1),dtype=float)
    
    for i in range(net.no_bus):
        net.Voltage2[i,0]=net.bus[i].Vm
        net.delta2[i,0]=net.bus[i].delta
        
    for i in range(net.no_fd):
        net.fd_P2[i,0]=net.fd[i].P[0]
        

# store_voltage2(net)



def calculate_aij(net):
    
    net.aij=np.zeros((net.no_fd,net.no_bus),dtype=float)
    
    X=net.Bd_inv

    for i in range(net.no_fd):
        
        n=net.fd[i].end_bus[0]
        m=net.fd[i].end_bus[1]
        
        for k in range(net.no_bus):
            
            net.aij[i,k]=(X[n,k]-X[m,k])/(net.fd[i].x)
    
            
calculate_aij(net)     

def calculate_fd_P_aij(net):
    
    net.fd_P_ail=np.zeros((net.no_fd,1),dtype=float)
    net.Pi=np.zeros((net.no_bus,1),dtype=float)
    
    for i in range(net.no_bus):
        net.Pi[i,0]=net.bus[i].P
        
    net.fd_P_ail=net.aij@net.Pi
    
calculate_fd_P_aij(net)   

def Precontingency_VoltageAngle(net):
    
    for i in range(net.no_bus):
        net.bus[i].Vm_pre=net.bus[i].Vm
        net.bus[i].delta_pre=net.bus[i].delta

#Precontingency_VoltageAngle(net)
        
    
    