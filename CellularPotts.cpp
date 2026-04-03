#include "CellularPotts.hpp"

template<typename T>
Array2d<T>::Array2d() : width{0}, height{0},data{0}
{}

template<typename T>
Array2d<T>::Array2d(size_t width, size_t height) : width{width},height{height},data{ static_cast<short unsigned int>(width*height)}
{}

template<typename T>
Array2d<T>::Array2d(size_t width, size_t height, T init) : width{width},height{height},data{}
{
        data.resize(width*height,init);
}

template<typename T>
const T& Array2d<T>::operator()(size_t x, size_t y) const{
        return data[x*width+y];
}
template<typename T>
T& Array2d<T>::operator()(size_t x, size_t y){
        return data[x*width+y];
}


Cell::Cell(uint16_t cell_id) :
        cell_id{cell_id},
        invasive{false},
        area{0},
        perim{0},
        perim2{0},
        adhesion_to_non_invasive{0},
        adhesion_to_invasive{0},
        last_point_x{0},
        last_point_y{0}
{}

CellularPotts::CellularPotts(size_t width, size_t height) : lattice{width,height,EMPTY}
{}

size_t CellularPotts::sample_int(size_t max) {
        return _sample_int(r_gen,std::uniform_int_distribution<size_t>::param_type {0,max});
}

size_t CellularPotts::sample_x(){
        return sample_int(lattice.width-1);
}

size_t CellularPotts::sample_y() {
        return sample_int(lattice.height-1);
}

bool CellularPotts::sample_bool(){
        return sample_uniform(r_gen) < 0.5;
}

uint16_t CellularPotts::sample_neighbor_state(size_t x, size_t y){
        if (sample_bool()){
                if (x == 0) return lattice(x+1,y);
                if (x == lattice.width-1) return lattice(x-1,y);
                if (sample_bool()) return lattice(x+1,y);
                return lattice(x-1,y);
        }
        else {
                if (y == 0) return lattice(x,y+1);
                if (y == lattice.height-1) return lattice(x,y-1);
                if (sample_bool()) return lattice(x,y+1);
                return lattice(x,y-1);
        }
}

void CellularPotts::remove_perim2(size_t x, size_t y){
        if (
                lattice(x,y) != EMPTY &&
                number_different_neighbor(x,y,lattice(x,y)) != 0
                ) cells[lattice(x,y)].perim2--;
        if (
                x!=0 && 
                lattice(x-1,y) != EMPTY &&
                number_different_neighbor(x-1,y,lattice(x-1,y)) != 0
                ) cells[lattice(x-1,y)].perim2--;
        if (
                y!=0 &&
                lattice(x,y-1) != EMPTY &&
                number_different_neighbor(x,y-1,lattice(x,y-1)) != 0
                ) cells[lattice(x,y-1)].perim2--;
        if (
                x!=lattice.width-1 &&
                lattice(x+1,y) != EMPTY &&
                number_different_neighbor(x+1,y,lattice(x+1,y)) != 0
                ) cells[lattice(x+1,y)].perim2--;
        if (
                y!=lattice.height-1 &&
                lattice(x,y+1) != EMPTY &&
                number_different_neighbor(x,y+1,lattice(x,y+1)) != 0
                ) cells[lattice(x,y+1)].perim2--;
}
void CellularPotts::add_perim2(size_t x, size_t y){
        if (
                lattice(x,y) != EMPTY &&
                number_different_neighbor(x,y,lattice(x,y)) != 0
                ) cells[lattice(x,y)].perim2++;
        if (
                x!=0 && 
                lattice(x-1,y) != EMPTY &&
                number_different_neighbor(x-1,y,lattice(x-1,y)) != 0
                ) cells[lattice(x-1,y)].perim2++;
        if (
                y!=0 &&
                lattice(x,y-1) != EMPTY &&
                number_different_neighbor(x,y-1,lattice(x,y-1)) != 0
                ) cells[lattice(x,y-1)].perim2++;
        if (
                x!=lattice.width-1 &&
                lattice(x+1,y) != EMPTY &&
                number_different_neighbor(x+1,y,lattice(x+1,y)) != 0
                ) cells[lattice(x+1,y)].perim2++;
        if (
                y!=lattice.height-1 &&
                lattice(x,y+1) != EMPTY &&
                number_different_neighbor(x,y+1,lattice(x,y+1)) != 0
                ) cells[lattice(x,y+1)].perim2++;
}

void CellularPotts::remove_adh(size_t x, size_t y){
        uint16_t cell_id = lattice(x,y);
        if (cell_id == EMPTY) return;
        uint16_t adhesion_to_non_invasive = 0;
        uint16_t adhesion_to_invasive = 0;
        number_adh(x,y,cell_id,adhesion_to_non_invasive,adhesion_to_invasive);
        cells[cell_id].adhesion_to_non_invasive -= adhesion_to_non_invasive;
        cells[cell_id].adhesion_to_invasive -= adhesion_to_invasive;
}

void CellularPotts::add_adh(size_t x, size_t y){
        uint16_t cell_id = lattice(x,y);
        if (cell_id == EMPTY) return;
        uint16_t adhesion_to_non_invasive = 0;
        uint16_t adhesion_to_invasive = 0;
        number_adh(x,y,cell_id,adhesion_to_non_invasive,adhesion_to_invasive);
        cells[cell_id].adhesion_to_non_invasive += adhesion_to_non_invasive;
        cells[cell_id].adhesion_to_invasive += adhesion_to_invasive;
}
void CellularPotts::remove_adh_zone(size_t x, size_t y){
        remove_adh(x,y);
        if (x!=0) remove_adh(x-1,y);
        if (y!=0) remove_adh(x,y-1);
        if (x!=lattice.width-1) remove_adh(x+1,y);
        if (x!=lattice.height-1) remove_adh(x,y+1);
}
void CellularPotts::add_adh_zone(size_t x, size_t y){
        add_adh(x,y);
        if (x!=0) add_adh(x-1,y);
        if (y!=0) add_adh(x,y-1);
        if (x!=lattice.width-1) add_adh(x+1,y);
        if (x!=lattice.height-1) add_adh(x,y+1);
}

void CellularPotts::update_lattice(size_t x, size_t y,uint16_t new_state){
        uint16_t old_state = lattice(x,y);
        // update area : one unit given from old_state to new_state
        if (old_state != EMPTY) cells[old_state].area--;
        if (new_state != EMPTY) cells[new_state].area++;
        // update perim : 4, 2, 0, -1 or -4 units given from old_to new_state
        uint16_t perim_point_old = number_different_neighbor(x,y,old_state);
        int16_t perim_given = 4 - (2*perim_point_old);
        if (old_state != EMPTY) cells[old_state].perim -= perim_given;
        if (new_state != EMPTY) cells[new_state].perim += perim_given;

        remove_perim2(x,y);
        remove_adh_zone(x,y);
        //update lattice
        lattice(x,y) = new_state;
        add_perim2(x,y);
        add_adh_zone(x,y);
}

void CellularPotts::MH_step(){
        size_t x = sample_x();
        size_t y = sample_y();
        uint16_t current_state = lattice(x,y);
        if (current_state == 0){
                // asm("int3");
        }
        uint16_t new_state = sample_neighbor_state(x,y);
        update_lattice(x,y,new_state);
        double H_new = compute_energy();
        double delta_H = H_new - H;
        double r = sample_uniform(r_gen);
        if (r >= std::exp(- delta_H/T)){
                // dont't accept transition : revert state
                update_lattice(x,y,current_state);
        }
        else {
                H = H_new;
        }
}

uint16_t CellularPotts::number_different_neighbor(size_t x, size_t y, uint16_t cell_id){
        uint16_t nb = 0;
        if (x == 0 || lattice(x-1,y) != cell_id) nb++;
        if (x == lattice.width-1 || lattice(x+1,y) != cell_id) nb++;
        if (y == 0 || lattice(x,y-1) != cell_id) nb++;
        if (y == lattice.height-1 || lattice(x,y+1) != cell_id) nb++;
        return nb;
}

void CellularPotts::number_adh(size_t x, size_t y, uint16_t cell_id, uint16_t& adhesion_to_non_invasive, uint16_t& adhesion_to_invasive){
        uint16_t nb = 0;
        if (x != 0 && lattice(x-1,y) != cell_id && lattice(x-1,y) != EMPTY){
                if (cells[lattice(x-1,y)].invasive){
                        adhesion_to_invasive++;
                }else {
                        adhesion_to_non_invasive++;
                }
        }
        if (x != lattice.width-1 && lattice(x+1,y) != cell_id && lattice(x+1,y) != EMPTY){
                if (cells[lattice(x+1,y)].invasive){
                        adhesion_to_invasive++;
                }else {
                        adhesion_to_non_invasive++;
                }
        }
        if (y != 0 && lattice(x,y-1) != cell_id && lattice(x,y-1) != EMPTY){
                if (cells[lattice(x,y-1)].invasive){
                        adhesion_to_invasive++;
                }else {
                        adhesion_to_non_invasive++;
                }
        }
        if (y != lattice.height && lattice(x,y+1) != cell_id && lattice(x,y+1) != EMPTY){
                if (cells[lattice(x,y+1)].invasive){
                        adhesion_to_invasive++;
                }else {
                        adhesion_to_non_invasive++;
                }
        }
}

void CellularPotts::initialize_cells_attributes(){
        for (Cell &cell : cells){
                initialize_cell_attributes(cell);
        }
}
void CellularPotts::initialize_cell_attributes(Cell &cell){
        cell.area = 0;
        cell.perim = 0;
        for(size_t x = 0; x<lattice.width; x++){
                for (size_t y = 0; y<lattice.height; y++){
                        if (lattice(x,y) == cell.cell_id){
                                cell.last_point_x = x;
                                cell.last_point_y = y;
                                cell.area++;
                                uint16_t nb_diff_neigh = number_different_neighbor(x,y,cell.cell_id);
                                cell.perim +=  nb_diff_neigh;
                                if (nb_diff_neigh != 0){
                                        cell.perim2 ++;
                                }
                                number_adh(x,y,cell.cell_id,cell.adhesion_to_non_invasive,cell.adhesion_to_invasive);
                                return;
                        }
                }
        }
}

double CellularPotts::compute_area_energy(){
        double H_area = 0;
        for (Cell& cell : cells){
                double area_diff = cell.area - target_area;
                H_area = area_diff*area_diff;
        }
        return lambda_area * H_area;
}

double CellularPotts::compute_perim_energy(){
        double H_perim = 0;
        for (Cell& cell : cells){
                double perim_diff = cell.perim2 - target_perim;
                H_perim = perim_diff*perim_diff;
        }
        return lambda_perim * H_perim;
}

double CellularPotts::compute_energy(){
        double H_area = compute_area_energy();
        double H_perim = compute_perim_energy();

        return H_area + H_perim;
}

