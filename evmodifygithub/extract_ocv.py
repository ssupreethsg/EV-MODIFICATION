import pybamm
import inspect

def main():
    param = pybamm.ParameterValues("Chen2020")
    
    print("Negative electrode OCP [V]:")
    print(inspect.getsource(param["Negative electrode OCP [V]"]))
    
    print("Positive electrode OCP [V]:")
    print(inspect.getsource(param["Positive electrode OCP [V]"]))
    
if __name__ == "__main__":
    main()
