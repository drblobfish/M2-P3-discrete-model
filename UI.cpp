#include <fstream>
#include <iostream>

#include "raylib.h"
#include "CellularPotts.hpp"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#define PADDING 20

int screenWidth = 800;
int screenHeight = 450;
float box_size = 0;
CellularPotts cp {100,100};
Rectangle cells_rec = {0};
Rectangle gui_rec = {0};

float iter_per_frame_f = 5000;
int iter_per_frame = 5000;
char iter_per_frame_str[10] = "";
bool running = true;
bool step_by_step = false;
bool reset = true;
int total_nb_iter = 0;

bool cell_visible = true;

uint16_t cell_id_focus = CellularPotts::EMPTY;
char cell_id_focus_info_str[100] = "";

char *param_file = nullptr;

RenderTexture2D cell_texture;

void draw_cells_texture(void){
        BeginTextureMode(cell_texture);
        for (size_t x = 0; x < cp.lattice.width; x++){
                for (size_t y = 0; y < cp.lattice.height; y++){
                        if (cell_visible && cp.lattice(x,y) != CellularPotts::EMPTY){
                                DrawRectangle(x,y,1,1,cp.cells[cp.lattice(x,y)].invasive ? ORANGE : PINK);
                        }
                        else {
                                DrawRectangle(x,y,1,1,
                                        ColorLerp(BLACK,SKYBLUE,5*cp.S_concentration(x,y)));
                        }
                }
        }
        EndTextureMode();
}
void draw_cell_walls(void){
        for (size_t x = 0; x < cp.lattice.width; x++){
                for (size_t y = 0; y < cp.lattice.height; y++){
                        if (x != cp.lattice.width-1 &&  cp.lattice(x,y) != cp.lattice(x+1,y)){
                                DrawLineEx(
                                        Vector2{.x = (x+1)*box_size, .y = y*box_size},
                                        Vector2{.x = (x+1)*box_size, .y = (y+1)*box_size},
                                        2, RED);
                        }
                        if (y != cp.lattice.height-1 && cp.lattice(x,y) != cp.lattice(x,y+1)){
                                DrawLineEx(
                                        Vector2{.x = (x)*box_size, .y = (y+1)*box_size},
                                        Vector2{.x = (x+1)*box_size, .y = (y+1)*box_size},
                                        2, RED);
                        }
                }
        }
}

void process_click_on_cell(void){
        if (!IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) return;
        Vector2 xy = GetMousePosition();
        if (!CheckCollisionPointRec(xy,cells_rec)) return;
        size_t x = (size_t) ((xy.x-cells_rec.x)/box_size);
        size_t y = (size_t) ((xy.y-cells_rec.y)/box_size);
        cell_id_focus = cp.lattice(x,y);
}

void draw_cell_info(void){
        if (cell_id_focus == CellularPotts::EMPTY) return;
        sprintf(cell_id_focus_info_str,"Cell %d\n%s\narea = %d\nperimeter = %d\nperimeter2 = %d\nadh_non_inv = %d\nadh_inv = %d\nS_concentration = %lf",
                        cell_id_focus,
                        cp.cells[cell_id_focus].invasive ? "invasive" : "non invasive",
                        cp.cells[cell_id_focus].area,
                        cp.cells[cell_id_focus].perim1,
                        cp.cells[cell_id_focus].perim2,
                        cp.cells[cell_id_focus].adhesion_to_non_invasive,
                        cp.cells[cell_id_focus].adhesion_to_invasive,
                        cp.cells[cell_id_focus].S_concentration
                        );
        GuiLabel((Rectangle){ gui_rec.x, gui_rec.y+4*PADDING, gui_rec.width, 5*PADDING }, cell_id_focus_info_str);
}

void parse_param_file(){
        if (param_file == nullptr) return;
        std::ifstream f(param_file);
        std::string param;
        while (f >> param){
                if (param == "T") f >> cp.T;
                if (param == "target_area") f >> cp.target_area;
                if (param == "lambda_area") f >> cp.lambda_area;
                if (param == "target_perim") f >> cp.target_perim;
                if (param == "lambda_perim") f >> cp.lambda_perim;
                if (param == "J_ii") f >> cp.J_ii;
                if (param == "J_in") f >> cp.J_in;
                if (param == "J_nn") f >> cp.J_nn;
                if (param == "S_diffusion") f >> cp.S_diffusion;
                if (param == "S_decay") f >> cp.S_decay;
                if (param == "S_emit") f >> cp.S_emit;
                if (param == "t_ni") f >> cp.t_ni;
                if (param == "t_in") f >> cp.t_in;
                if (param == "mu_S_inv") f >> cp.mu_S_inv;
                if (param == "xi") f >> cp.xi;
                if (param == "init_with_I") f >> cp.init_with_I;
                if (param == "init_R") f >> cp.init_R;
                if (param == "init_nb_cells") f >> cp.init_nb_cells;
        }
}

int main(int argc, char*argv[])
{
        if (argc == 2) param_file = argv[1];
        parse_param_file();
        InitWindow(screenWidth, screenHeight, "Cellular Potts Model");
        SetTargetFPS(60);
        cell_texture = LoadRenderTexture(cp.lattice.width, cp.lattice.height);

        while (!WindowShouldClose())
        {
                if (reset){
                        cp.initialize_board();
                        reset = false;
                }
                if (running){
                        for (size_t i = 0; i<iter_per_frame;i++){
                                cp.MH_step();
                        }
                        if (step_by_step) {
                                running = false;
                                step_by_step = false;
                        }
                        total_nb_iter += iter_per_frame;
                }

                screenWidth = GetScreenWidth();
                screenHeight = GetScreenHeight();
                float size_x = (float) screenWidth / cp.lattice.width;
                float size_y = (float) screenHeight / cp.lattice.height;
                box_size = (size_x < size_y) ? size_x : size_y;
                cells_rec = Rectangle{.x = 0, .y = 0, .width = cp.lattice.width * box_size, .height = cp.lattice.height*box_size};
                gui_rec = Rectangle{.x = cells_rec.width+PADDING, .y = PADDING, .width = screenWidth - cells_rec.width - 2 * PADDING, .height = (float) screenHeight - PADDING};

                // Draw
                draw_cells_texture();

                BeginDrawing();
                ClearBackground(RAYWHITE);
                process_click_on_cell();
                DrawRectangleLinesEx(cells_rec, 2.0, BLACK);
                reset = GuiButton((Rectangle){ gui_rec.x, gui_rec.y, 3*PADDING, PADDING }, "reset");
                if (running){
                        running = !GuiButton((Rectangle){ gui_rec.x+4*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_PLAYER_PAUSE,nullptr));
                }else{
                        running = GuiButton((Rectangle){ gui_rec.x+4*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_PLAYER_PLAY,nullptr));
                }
                if(GuiButton((Rectangle){ gui_rec.x + 6*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_PLAYER_NEXT,nullptr))){
                        step_by_step = true;
                        running = true;
                }
                GuiToggle((Rectangle){ gui_rec.x + 8*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_EYE_ON,nullptr),&cell_visible);
                GuiSlider(Rectangle{gui_rec.x+ 60,gui_rec.y+ 2* PADDING,gui_rec.width-80,PADDING}, "iter/frame",iter_per_frame_str, &iter_per_frame_f, 1, 5e3);
                iter_per_frame = (int) iter_per_frame_f;
                sprintf(iter_per_frame_str,"%d",iter_per_frame);
                draw_cell_info();
                DrawTexturePro(
                        cell_texture.texture,
                        Rectangle{0,0,(float)cell_texture.texture.width,-(float)cell_texture.texture.height},
                        cells_rec,
                        Vector2{0,0},
                        0,
                        WHITE);
                draw_cell_walls();
                std::cout << "iter = "<< total_nb_iter <<"\n";

                EndDrawing();
        }
        UnloadRenderTexture(cell_texture);
        CloseWindow();

        return 0;
}
