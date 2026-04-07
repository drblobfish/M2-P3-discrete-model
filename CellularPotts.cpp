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
        perim1{0},
        perim2{0},
        adhesion_to_non_invasive{0},
        adhesion_to_invasive{0},
        S_concentration{0},
        last_point_x{0},
        last_point_y{0}
{}

CellularPotts::CellularPotts(size_t width, size_t height) : lattice{width,height,EMPTY},S_concentration{width,height,0},S_concentration_back{width,height,0}
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

uint16_t CellularPotts::sample_neighbor_state(size_t x, size_t y, size_t& nx, size_t& ny){
        if (sample_bool()){
                if (x == 0) {nx = x+1; ny = y;}
                else if (x == lattice.width-1) {nx = x-1; ny = y;}
                else if (sample_bool()) {nx = x+1; ny = y;}
                else {nx = x-1; ny = y;}
        }
        else {
                if (y == 0) {nx = x; ny = y+1;}
                else if (y == lattice.height-1) {nx = x; ny = y-1;}
                else if (sample_bool()) {nx = x; ny = y+1;}
                else {nx = x; ny = y-1;}
        }
        return lattice(nx,ny);
}

void CellularPotts::initialize_board(){
        initialize_S_concentration();
        initialize_tumor_core(lattice.width/2,lattice.height/2,10,20);
        // lattice(50,50) = 0;
        // lattice(60,60) = 1;
        // cells.push_back(Cell(0));
        // cells.push_back(Cell(1));
}

void CellularPotts::initialize_S_concentration(){
        for (size_t x = 0; x<lattice.width;x++){
                for (size_t y = 0; y<lattice.height; y++){
                        S_concentration(x,y) = 0;
                }
        }
}


void CellularPotts:: initialize_tumor_core(double center_x, double center_y, double radius, size_t nb_cells){
        for (size_t x = 0; x<lattice.width;x++){
                for (size_t y = 0; y<lattice.height; y++){
                        double dx = center_x - x;
                        double dy = center_y - y;
                        double r = std::sqrt(dx*dx + dy*dy);
                        lattice(x,y) = r<=radius ? 0 : EMPTY;
                }
        }
        std::vector<double> meanx = std::vector<double>(nb_cells);
        std::vector<double> meany = std::vector<double>(nb_cells);
        for (size_t i = 0; i<nb_cells; i++){
                do {
                        meanx[i] = (sample_uniform(r_gen)-0.5)*2*radius + center_x;
                        meany[i] = (sample_uniform(r_gen)-0.5)*2*radius + center_y;
                }while (std::sqrt((center_x - meanx[i])*(center_x - meanx[i])+(center_y - meany[i])*(center_y - meany[i])) > radius);
        }
        for (size_t i = 0; i< 100; i++){
                for (size_t x = 0; x<lattice.width;x++){
                        for (size_t y = 0; y<lattice.height; y++){
                                if (lattice(x,y) != EMPTY){
                                        uint16_t best_cell_id = 0;
                                        double best_dist = 2*radius;
                                        for (uint16_t j = 0; j<nb_cells; j++){
                                                double dist = std::sqrt((x-meanx[j])*(x-meanx[j]) + (y-meany[j])*(y-meany[j]));
                                                if (dist < best_dist){
                                                        best_cell_id = j;
                                                        best_dist = dist;
                                                }
                                        }
                                        lattice(x,y) = best_cell_id;
                                }
                        }
                }
                for (size_t j=0;j<nb_cells;j++){
                        size_t N = 0;
                        meanx[j] = 0;
                        meany[j] = 0;
                        for (size_t x = 0; x<lattice.width;x++){
                                for (size_t y = 0; y<lattice.height; y++){
                                        if (lattice(x,y) == j){
                                                N++;
                                                meanx[j] += x;
                                                meany[j] += y;
                                        }
                                }
                        }
                        meanx[j] /= N;
                        meany[j] /= N;
                }
        }
        cells.clear();
        for (uint16_t i=0; i<nb_cells; i++){
                cells.push_back(Cell(i));
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

void CellularPotts::remove_adh(size_t x1, size_t y1,size_t x2,size_t y2){
        uint16_t cell_id1 = lattice(x1,y1);
        uint16_t cell_id2 = lattice(x2,y2);
        if (cell_id1 == cell_id2) return;
        if (cell_id1 != EMPTY) cells[cell_id1].perim1--;
        if (cell_id2 != EMPTY) cells[cell_id2].perim1--;
        if (cell_id1 == EMPTY || cell_id2 == EMPTY) return;
        if (cells[cell_id1].invasive){
                cells[cell_id2].adhesion_to_invasive--;
        }else {
                cells[cell_id2].adhesion_to_non_invasive--;
        }
        if (cells[cell_id2].invasive){
                cells[cell_id1].adhesion_to_invasive--;
        }else {
                cells[cell_id1].adhesion_to_non_invasive--;
        }
}

void CellularPotts::add_adh(size_t x1, size_t y1,size_t x2,size_t y2){
        uint16_t cell_id1 = lattice(x1,y1);
        uint16_t cell_id2 = lattice(x2,y2);
        if (cell_id1 == cell_id2) return;
        if (cell_id1 != EMPTY) cells[cell_id1].perim1++;
        if (cell_id2 != EMPTY) cells[cell_id2].perim1++;
        if (cell_id1 == EMPTY || cell_id2 == EMPTY) return;
        if (cells[cell_id1].invasive){
                cells[cell_id2].adhesion_to_invasive++;
        }else {
                cells[cell_id2].adhesion_to_non_invasive++;
        }
        if (cells[cell_id2].invasive){
                cells[cell_id1].adhesion_to_invasive++;
        }else {
                cells[cell_id1].adhesion_to_non_invasive++;
        }
}

void CellularPotts::remove_adh_zone(size_t x, size_t y){
        if (x!=0) remove_adh(x,y,x-1,y);
        if (y!=0) remove_adh(x,y,x,y-1);
        if (x!=lattice.width-1) remove_adh(x,y,x+1,y);
        if (x!=lattice.height-1) remove_adh(x,y,x,y+1);
}
void CellularPotts::add_adh_zone(size_t x, size_t y){
        if (x!=0) add_adh(x,y,x-1,y);
        if (y!=0) add_adh(x,y,x,y-1);
        if (x!=lattice.width-1) add_adh(x,y,x+1,y);
        if (x!=lattice.height-1) add_adh(x,y,x,y+1);
}

void CellularPotts::update_lattice(size_t x, size_t y,uint16_t new_state){
        uint16_t old_state = lattice(x,y);
        // update area : one unit given from old_state to new_state
        if (old_state != EMPTY) cells[old_state].area--;
        if (new_state != EMPTY) cells[new_state].area++;
        // update perim : 4, 2, 0, -1 or -4 units given from old_to new_state
        // uint16_t perim_point_old = number_different_neighbor(x,y,old_state);
        // int16_t perim_given = 4 - (2*perim_point_old);
        // if (old_state != EMPTY) cells[old_state].perim -= perim_given;
        // if (new_state != EMPTY) cells[new_state].perim += perim_given;

        remove_perim2(x,y);
        remove_adh_zone(x,y);
        remove_S_concentration(x,y);
        //update lattice
        lattice(x,y) = new_state;
        add_perim2(x,y);
        add_adh_zone(x,y);
        add_S_concentration(x,y);
}

void CellularPotts::MH_step(){
        size_t x = sample_x();
        size_t y = sample_y();
        uint16_t current_state = lattice(x,y);
        // asm("int3");
        size_t neighbor_x;
        size_t neighbor_y;
        uint16_t new_state = sample_neighbor_state(x,y,neighbor_x,neighbor_y);
        update_lattice(x,y,new_state);
        double H_new = compute_energy();
        double chemotactic_bias = xi * (S_concentration(x,y) - S_concentration(neighbor_x,neighbor_y));
        double delta_H = H_new - H;
        if (new_state != EMPTY && cells[new_state].invasive) delta_H += chemotactic_bias;
        double r = sample_uniform(r_gen);
        if (r >= std::exp(- delta_H/T)){
                // dont't accept transition : revert state
                update_lattice(x,y,current_state);
        }
        else {
                H = H_new;
        }
        if ((S_step_counter++) == S_step_freq){
                S_step_counter = 0;
                S_step();
                update_cell_S_concetration();
        }
        CTMC_step();
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
        cell.perim1 = 0;
        cell.perim2 = 0;
        cell.adhesion_to_non_invasive = 0;
        cell.adhesion_to_invasive = 0;
        for(size_t x = 0; x<lattice.width; x++){
                for (size_t y = 0; y<lattice.height; y++){
                        if (lattice(x,y) == cell.cell_id){
                                cell.last_point_x = x;
                                cell.last_point_y = y;
                                cell.area++;
                                uint16_t nb_diff_neigh = number_different_neighbor(x,y,cell.cell_id);
                                cell.perim1 +=  nb_diff_neigh;
                                if (nb_diff_neigh != 0){
                                        cell.perim2 ++;
                                }
                                number_adh(x,y,cell.cell_id,cell.adhesion_to_non_invasive,cell.adhesion_to_invasive);
                        }
                }
        }
}

double CellularPotts::compute_area_energy(){
        double H_area = 0;
        for (Cell& cell : cells){
                double area_diff = cell.area - target_area;
                H_area += area_diff*area_diff;
        }
        return lambda_area * H_area;
}

double CellularPotts::compute_perim_energy(){
        double H_perim = 0;
        for (Cell& cell : cells){
                double perim_diff = cell.perim2 - target_perim;
                H_perim += perim_diff*perim_diff;
        }
        return lambda_perim * H_perim;
}

double CellularPotts::compute_adhesion_energy(){
        double H_adh = 0;
        for (Cell& cell : cells){
                if (cell.invasive){
                        H_adh += J_ii * cell.adhesion_to_invasive + J_in * cell.adhesion_to_non_invasive;
                } else {
                        H_adh += J_in * cell.adhesion_to_invasive + J_nn * cell.adhesion_to_non_invasive;
                }
        }
        return H_adh;
}

double CellularPotts::compute_chemo_energy(){
        double H_chemo = 0;
        for (Cell& cell : cells){
                if (cell.invasive){
                        H_chemo += mu_S_inv * cell.S_concentration / cell.area;
                }
        }
        return H_chemo;
}

double CellularPotts::compute_energy(){
        double H_area = compute_area_energy();
        double H_perim = compute_perim_energy();
        double H_adh = compute_adhesion_energy();
        double H_chemo = compute_chemo_energy();

        return H_area + H_perim + H_adh + H_chemo;
}

void CellularPotts::S_step_point(size_t x, size_t y){
        double change = 0;
        if (x !=0) change += S_concentration(x-1,y);
        if (y !=0) change += S_concentration(x,y-1);
        if (x !=lattice.width-1) change += S_concentration(x+1,y);
        if (y !=lattice.height-1) change += S_concentration(x,y+1);
        change += - 4*S_concentration(x,y);
        change *= S_diffusion;
        change += - S_decay * S_concentration(x,y);
        if (lattice(x,y) != EMPTY && (!cells[lattice(x,y)].invasive)) change += S_emit;
        S_concentration_back(x,y) = S_concentration(x,y) + S_dt * change;
        if (S_concentration_back(x,y) > S_concentration_max) S_concentration_max = S_concentration_back(x,y);
}

/*
void CellularPotts::S_step_point(size_t x, size_t y){
        double dx = x - 0.5 * lattice.width;
        double dy = y - 0.5 * lattice.height;
        double d = std::sqrt(dx*dx + dy*dy);
        S_concentration_back(x,y) = 0.2 - 0.005 * d;
}
*/

void CellularPotts::S_step(){
        S_concentration_max = 0;
        for(size_t x = 0; x<lattice.width; x++){
                for (size_t y = 0; y<lattice.height; y++){
                        S_step_point(x,y);
                }
        }
        printf("max S concentration = %f\n",S_concentration_max);
        std::swap(S_concentration.data,S_concentration_back.data);
}

void CellularPotts::CTMC_step_cell(Cell& cell){
        double r = sample_uniform(r_gen);
        double connexion_ratio = (cell.adhesion_to_non_invasive + cell.adhesion_to_invasive)/(double) cell.perim1;
        if (cell.invasive){
                if (r<t_in * CTMC_dt * connexion_ratio){
                        cell.invasive = false;
                        initialize_cells_attributes();
                }
        } else {
                if (r<t_ni * CTMC_dt * (1-connexion_ratio)){
                        cell.invasive = true;
                        initialize_cells_attributes();
                }
        }
}

void CellularPotts::CTMC_step(){
        for (Cell& cell : cells){
                CTMC_step_cell(cell);
        }
}

void CellularPotts::update_cell_S_concetration(){
        for (Cell& cell : cells) cell.S_concentration = 0;
        for(size_t x = 0; x<lattice.width; x++){
                for (size_t y = 0; y<lattice.height; y++){
                        if (lattice(x,y) != EMPTY){
                                cells[lattice(x,y)].S_concentration += S_concentration(x,y);
                        }
                }
        }
}
void CellularPotts::remove_S_concentration(size_t x, size_t y){
        if (lattice(x,y) == EMPTY) return;
        cells[lattice(x,y)].S_concentration -= S_concentration(x,y);
}
void CellularPotts::add_S_concentration(size_t x, size_t y){
        if (lattice(x,y) == EMPTY) return;
        cells[lattice(x,y)].S_concentration += S_concentration(x,y);
}
