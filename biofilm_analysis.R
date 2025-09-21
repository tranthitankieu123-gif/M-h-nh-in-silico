
library(BacArena)
library(ggplot2)
library(dplyr)
library(viridis)

# Đọc dữ liệu mô phỏng
load("biomass_check_1_2_sim.RData")
# load("sim_oxy_v3_chemotaxis.RData")
# load("sim_after_12h_v2_20250920_122718.RData")

simulation_result <- sim2;
# str(simulation_result@specs)
# print(slotNames(simulation_result@specs))

# ===================================================================
pdf("sim_after_12h_v2_20250920_122718_analysis.pdf")
biofilm_data = evalArena(simulation_result, legend_pos = "bottom", retdata= TRUE, time=length(simulation_result@medlist)-1)
plotSubUsage(simulation_result, ret_data = FALSE) # Oxy
plotCurves(simulation_result, medplot=c("EX_cpd00007_e0"))
subUsage = plotSubUsage(simulation_result, ret_data = TRUE) # Oxy
plotSubDist2(simulation_result, sub="EX_cpd00007_e0") # Oxy
print(subUsage)
grid_size <- simulation_result@n
lastDataFrame <- biofilm_data$Population$time20


# Tạo một data frame chứa số lượng tế bào tại mỗi tọa độ (x, y)
density_total <- lastDataFrame %>%
  group_by(x, y) %>%
  summarise(TotalCell = n(), .groups = 'drop')

# Vẽ bản đồ nhiệt
heatmap_total <- ggplot(density_total, aes(x = x, y = y, fill = TotalCell)) +
  geom_hex() +
  coord_fixed() +
  scale_fill_binned(
    name = "Tong Sinh khoi\n(nM)",
    type = "viridis", # Sử dụng bảng màu viridis (xanh-vàng)
    option = "D",     # Tùy chọn "D" (viridis) là thang màu từ xanh dương-xanh lá-vàng
    n.breaks = 7,     # Chia dải màu thành 7 bậc để phân biệt rõ
    na.value = "grey50",
    guide = guide_coloursteps(even.steps = TRUE, show.limits = TRUE) # Hiển thị các bậc màu rõ ràng trên chú giải
  )+
  scale_fill_viridis() + 
  theme_minimal() +
  labs(
    title = "Biofilm by Total Cell Density",
    x = "X",
    y = "Y"
  )

print(heatmap_total)


biomass_density <- lastDataFrame %>%
  group_by(x, y) %>%
  summarise(TotalBiomass = sum(biomass), .groups = 'drop')

heatmap_biomass <- ggplot(biomass_density, aes(x = x, y = y, fill = TotalBiomass)) +
  geom_tile() +
  scale_fill_binned(
    name = "Tong Sinh khoi\n(nM)",
    type = "viridis", # Sử dụng bảng màu viridis (xanh-vàng)
    option = "D",     # Tùy chọn "D" (viridis) là thang màu từ xanh dương-xanh lá-vàng
    n.breaks = 7,     # Chia dải màu thành 7 bậc để phân biệt rõ
    na.value = "grey50",
    guide = guide_coloursteps(even.steps = TRUE, show.limits = TRUE) # Hiển thị các bậc màu rõ ràng trên chú giải
  ) +
  coord_fixed() +
  labs(
    title = "Bản đồ nhiệt Phân bố Mật độ Sinh khối Biofilm",
    subtitle = "Màu sắc thể hiện tổng khối lượng vật chất sống tại mỗi vị trí",
    x = "Toa do X",
    y = "Toa do Y"
  ) 
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        text = element_text(family = "Arial"))

# Hiển thị biểu đồ
print(heatmap_biomass)


# Tạo data frame mật độ, nhóm theo 'type' (loài)
density_by_species <- lastDataFrame %>%
  group_by(type, x, y) %>%
  summarise(CellCount = n(), .groups = 'drop') %>%
  mutate(Species = factor(paste("Loài", type))) # Tạo cột mới cho tên loài dễ đọc hơn

# Vẽ bản đồ nhiệt, sử dụng facet_wrap để tách mỗi loài thành một biểu đồ con
heatmap_species <- ggplot(density_by_species, aes(x = x, y = y, fill = CellCount)) +
  geom_raster() +
  scale_fill_viridis(name = "Mật độ Tế bào", option = "C") +
  coord_fixed() +
  facet_wrap(~ Species) + # Đây là lệnh tạo ra các biểu đồ con theo loài
  labs(
    title = "Biofilm by Org",
    subtitle = "Cho thấy sự hình thành ổ sinh thái (niche) của từng loài",
    x = "Tọa độ X",
    y = "Tọa độ Y"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold")) # Chỉnh sửa tiêu đề của biểu đồ con

# Hiển thị biểu đồ
print(heatmap_species)


# ===================================================================
# BƯỚC 5: PHÂN TÍCH 3 - PHÂN BỐ CỦA CÁC KIỂU HÌNH TRAO ĐỔI CHẤT
# ===================================================================
# Mục đích: Đây là phân tích chức năng. Nó cho thấy "ai đang làm gì và ở đâu".
#           Ví dụ: các tế bào lên men có tập trung ở lõi thiếu oxy không?

# Tạo data frame mật độ, nhóm theo 'phenotype'
density_by_phenotype <- lastDataFrame %>%
  group_by(phenotype, x, y) %>%
  summarise(CellCount = n(), .groups = 'drop') %>%
  mutate(Phenotype = factor(paste("Phenon ", phenotype)))

# Vẽ bản đồ nhiệt, facet theo 'phenotype'
heatmap_phenotype <- ggplot(density_by_phenotype, aes(x = x, y = y, fill = CellCount)) +
  geom_tile() +
  scale_fill_viridis(name = "Mật độ Tế bào", option = "C") +
  coord_fixed() +
  facet_wrap(~ Phenotype) +
  labs(
    title = "Phân bố Không gian theo Kiểu hình Trao đổi chất (t = 12)",
    subtitle = "Tiết lộ cấu trúc chức năng của cộng đồng biofilm",
    x = "Tọa độ X",
    y = "Tọa độ Y"
  ) +
  # theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

# Hiển thị biểu đồ
print(heatmap_phenotype)