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

struct Cell{
        uint16_t type;
        uint16_t area;
        uint16_t perim;
        uint16_t perim2;
        size_t last_point_x;
        size_t last_point_y;

        Cell(uint16_t type);
};

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

        double target_area = 25;
        double lambda_area = 1;
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
        void remove_perim2(size_t x, size_t y);
        void add_perim2(size_t x, size_t y);
        void update_lattice(size_t x, size_t y,uint16_t new_state);
        uint16_t number_different_neighbor(size_t x, size_t y, uint16_t type);

};

