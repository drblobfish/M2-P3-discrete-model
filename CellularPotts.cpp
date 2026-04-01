#include <cstdint>
#include <vector>
#include <random>

template<typename T>
struct Array2d{
        size_t width;
        size_t height;
        std::vector<T> data;

        Array2d();
        Array2d(size_t width,size_t height);
        Array2d(size_t width,size_t height,T init);

        const T& operator()(size_t x, size_t y) const;
        T& operator()(size_t x, size_t y);
};

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

struct Cell{
        uint16_t type;
        uint16_t area;
        uint16_t perim;
        size_t last_point_x;
        size_t last_point_y;

        Cell(uint16_t type);
};

Cell::Cell(uint16_t type) : type{type}, area{0},perim{0},last_point_x{0},last_point_y{0}{
}


struct CellularPotts{
        Array2d<uint16_t> lattice;
 
        std::random_device random_dev;
        std::default_random_engine r_gen {random_dev()};
        std::uniform_real_distribution<double> sample_uniform {};
        std::uniform_int_distribution<size_t> _sample_int;

        double H;
        double T = 1;

        constexpr static uint16_t EMPTY = std::numeric_limits<uint16_t>::max();
        std::vector<Cell> cells;

        double target_area = 10;
        double lambda_area = 10;
        double target_perim = 10;
        double lambda_perim = 1;


        CellularPotts(size_t width, size_t height);
        size_t sample_int(size_t max);
        size_t sample_x();
        size_t sample_y();
        bool sample_bool();
        uint16_t sample_neighbor_state(size_t x, size_t y);
        void MH_step();

        void initialize_cells_attributes();
        void initialize_cell_attributes(Cell &cell);
        double compute_energy();
        double compute_area_energy();
        double compute_perim_energy();
        void update_lattice(size_t x, size_t y,uint16_t new_state);
        uint16_t number_different_neighbor(size_t x, size_t y, uint16_t type);

};

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
        lattice(x,y) = new_state;
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

uint16_t CellularPotts::number_different_neighbor(size_t x, size_t y, uint16_t type){
        uint16_t nb = 0;
        if (x == 0 || lattice(x-1,y) != type) nb++;
        if (x == lattice.width-1 || lattice(x+1,y) != type) nb++;
        if (y == 0 || lattice(x,y-1) != type) nb++;
        if (y == lattice.height-1 || lattice(x,y+1) != type) nb++;
        return nb;
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
                        if (lattice(x,y) == cell.type){
                                cell.last_point_x = x;
                                cell.last_point_y = y;
                                cell.area++;
                                cell.perim += number_different_neighbor(x,y,cell.type);
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
                double perim_diff = cell.perim - target_perim;
                H_perim = perim_diff*perim_diff;
        }
        return lambda_perim * H_perim;
}

double CellularPotts::compute_energy(){
        double H_area = compute_area_energy();
        double H_perim = compute_area_energy();

        return H_area + H_perim;
}

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
