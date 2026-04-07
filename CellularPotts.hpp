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
        uint16_t cell_id;
        bool invasive;
        uint16_t area;
        uint16_t perim1;
        uint16_t perim2;
        uint16_t adhesion_to_non_invasive;
        uint16_t adhesion_to_invasive;
        double S_concentration;
        double S_min_concentration;
        size_t last_point_x;
        size_t last_point_y;

        Cell(uint16_t cell_id);
};

struct CellularPotts{
        Array2d<uint16_t> lattice;
 
        std::random_device random_dev;
        std::default_random_engine r_gen {random_dev()};
        std::uniform_real_distribution<double> sample_uniform {};
        std::uniform_int_distribution<size_t> _sample_int;

        double H;
        double T = 5;

        constexpr static uint16_t EMPTY = std::numeric_limits<uint16_t>::max();
        std::vector<Cell> cells;

        double target_area = 25;
        double lambda_area = 1;
        double target_perim = 10;
        double lambda_perim = 1;
        double J_ii = 2;
        double J_in = 2;
        double J_nn = 1;

        Array2d<double> S_concentration;
        Array2d<double> S_concentration_back;
        double S_diffusion = 500;
        double S_decay = 1;
        double S_emit = 2;
        double S_dt = 5e-5;
        double S_concentration_max = 1;
        size_t S_step_freq = 100;
        size_t S_step_counter = 0;

        double t_ni = 0.1;
        double t_in = 0.01;
        double CTMC_dt = 1e-4;

        double mu_S_inv = 0;
        double xi = 40;

        CellularPotts(size_t width, size_t height);
        size_t sample_int(size_t max);
        size_t sample_x();
        size_t sample_y();
        bool sample_bool();
        uint16_t sample_neighbor_state(size_t x, size_t y, size_t& nx, size_t& ny);
        void MH_step();

        void initialize_board();
        void initialize_S_concentration();
        void initialize_tumor_core(double center_x, double center_y, double radius, size_t nb_cells);
        void initialize_cells_attributes();
        void initialize_cell_attributes(Cell &cell);
        double compute_energy();
        double compute_area_energy();
        double compute_perim_energy();
        double compute_adhesion_energy();
        double compute_chemo_energy();
        void remove_perim2(size_t x, size_t y);
        void add_perim2(size_t x, size_t y);
        void number_adh(size_t x, size_t y, uint16_t cell_id, uint16_t& adhesion_to_non_invasive, uint16_t& adhesion_to_invasive);
        void remove_adh(size_t x1, size_t y1,size_t x2, size_t y2);
        void remove_adh_zone(size_t x, size_t y);
        void add_adh(size_t x1, size_t y1,size_t x2, size_t y2);
        void add_adh_zone(size_t x, size_t y);
        void update_lattice(size_t x, size_t y,uint16_t new_state);
        uint16_t number_different_neighbor(size_t x, size_t y, uint16_t cell_id);

        void S_step();
        void S_step_point(size_t x, size_t y);

        void CTMC_step_cell(Cell& cell);
        void CTMC_step();
        void update_cell_S_concetration();
        void remove_S_concentration(size_t x, size_t y);
        void add_S_concentration(size_t x,size_t y);
};

