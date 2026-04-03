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

float iter_per_frame_f = 100;
int iter_per_frame = 100;
char iter_per_frame_str[10] = "";
bool running = false;
bool step_by_step = false;

uint16_t cell_id_focus = CellularPotts::EMPTY;
char cell_id_focus_info_str[100] = "";

void draw_cells(void){
        for (size_t x = 0; x < cp.lattice.width; x++){
                for (size_t y = 0; y < cp.lattice.height; y++){
                        if (cp.lattice(x,y) != CellularPotts::EMPTY)
                                DrawRectangle(x*box_size,y*box_size,box_size,box_size,PINK);
                }
        }
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
        sprintf(cell_id_focus_info_str,"Cell %d\narea = %d\nperimeter = %d\nadh_non_inv = %d\nadh_inv = %d",
                        cell_id_focus,
                        cp.cells[cell_id_focus].area,
                        cp.cells[cell_id_focus].perim2,
                        cp.cells[cell_id_focus].adhesion_to_non_invasive,
                        cp.cells[cell_id_focus].adhesion_to_invasive
                        );
        GuiLabel((Rectangle){ gui_rec.x, gui_rec.y+4*PADDING, gui_rec.width, 5*PADDING }, cell_id_focus_info_str);
}

int main(void)
{
        InitWindow(screenWidth, screenHeight, "Cellular Potts Model");
        SetTargetFPS(60);

        cp.initialize_board(cp.lattice.width/2,cp.lattice.height/2,10,20);
        cp.initialize_cells_attributes();
        cp.H = cp.compute_energy();

        while (!WindowShouldClose())
        {
                if (running){
                        for (size_t i = 0; i<iter_per_frame;i++){
                                cp.MH_step();
                        }
                        if (step_by_step) {
                                running = false;
                                step_by_step = false;
                        }
                }

                screenWidth = GetScreenWidth();
                screenHeight = GetScreenHeight();
                float size_x = (float) screenWidth / cp.lattice.width;
                float size_y = (float) screenHeight / cp.lattice.height;
                box_size = (size_x < size_y) ? size_x : size_y;
                cells_rec = Rectangle{.x = 0, .y = 0, .width = cp.lattice.width * box_size, .height = cp.lattice.height*box_size};
                gui_rec = Rectangle{.x = cells_rec.width+PADDING, .y = PADDING, .width = screenWidth - cells_rec.width - 2 * PADDING, .height = (float) screenHeight - PADDING};

                // Draw
                BeginDrawing();
                ClearBackground(RAYWHITE);
                process_click_on_cell();
                DrawRectangleLinesEx(cells_rec, 2.0, BLACK);
                GuiButton((Rectangle){ gui_rec.x, gui_rec.y, 3*PADDING, PADDING }, "reset");
                if (running){
                        running = !GuiButton((Rectangle){ gui_rec.x+4*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_PLAYER_PAUSE,nullptr));
                }else{
                        running = GuiButton((Rectangle){ gui_rec.x+4*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_PLAYER_PLAY,nullptr));
                }
                if(GuiButton((Rectangle){ gui_rec.x + 6*PADDING, gui_rec.y, PADDING, PADDING }, GuiIconText(ICON_PLAYER_NEXT,nullptr))){
                        step_by_step = true;
                        running = true;
                }
                GuiSlider(Rectangle{gui_rec.x+ 60,gui_rec.y+ 2* PADDING,gui_rec.width-60,PADDING}, "iter/frame",iter_per_frame_str, &iter_per_frame_f, 1, 1e4);
                iter_per_frame = (int) iter_per_frame_f;
                sprintf(iter_per_frame_str,"%d",iter_per_frame);
                draw_cell_info();
                draw_cells();
                draw_cell_walls();

                EndDrawing();
        }

        CloseWindow();

        return 0;
}
