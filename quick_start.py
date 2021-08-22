from common_modification_search import open_search 
from blind_search import blind_search 
from close_identify import close_search 
import os 


if __name__ == "__main__":
    current_path = os.getcwd() 
    print('Welcome to use pChem!')
    
    open_search(current_path)
    blind_search(current_path)
    close_search(current_path) 
    
    