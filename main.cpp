#include "CellularPotts.hpp"

int main(){
        CellularPotts cp {10,10};
        cp.lattice(5,5) = 0;
        cp.cells.push_back(Cell(0));
        cp.initialize_cells_attributes();
        cp.H = cp.compute_energy();
        for (size_t i = 0; i< 1000; i++){
                cp.MH_step();
        }
        return 1;
}
