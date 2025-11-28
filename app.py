import streamlit as st
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, LineString
from pyproj import Transformer
from scipy.interpolate import splprep, splev, interp1d
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from scipy.special import comb
import folium
from streamlit_folium import st_folium
from reportlab.lib.pagesizes import A4, letter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from io import BytesIO
from datetime import datetime

# 设置页面配置
st.set_page_config(page_title="黄岩岛领海基线面积计算", layout="wide")

st.title("黄岩岛领海基线面积计算")
st.markdown("基于2024年11月10日发布的《中华人民共和国政府关于黄岩岛领海基线的声明》")

# 定义领海基线标志点数据
baseline_data = [
    ("黄岩岛1", "15°08.1′", "117°50.9′"),
    ("黄岩岛2", "15°07.4′", "117°50.8′"),
    ("黄岩岛3", "15°07.0′", "117°50.6′"),
    ("黄岩岛4", "15°06.6′", "117°50.2′"),
    ("黄岩岛5", "15°06.1′", "117°49.5′"),
    ("黄岩岛6", "15°06.3′", "117°44.2′"),
    ("黄岩岛7", "15°07.3′", "117°43.1′"),
    ("黄岩岛8", "15°12.7′", "117°42.6′"),
    ("黄岩岛9", "15°13.1′", "117°42.8′"),
    ("黄岩岛10", "15°13.4′", "117°43.3′"),
    ("黄岩岛11", "15°13.5′", "117°43.9′"),
    ("黄岩岛12", "15°13.5′", "117°44.4′"),
    ("黄岩岛13", "15°09.6′", "117°49.7′"),
    ("黄岩岛14", "15°09.0′", "117°50.4′"),
    ("黄岩岛15", "15°08.5′", "117°50.8′"),
    ("黄岩岛1", "15°08.1′", "117°50.9′"),
]

def parse_coordinate(coord_str):
    """
    将度分格式（如 15°08.1′）或十进制度数转换为十进制度数
    """
    # 如果已经是数字类型，直接返回
    if isinstance(coord_str, (int, float)):
        return float(coord_str)
    
    # 转换为字符串
    coord_str = str(coord_str).strip()
    
    # 如果不包含度分符号，尝试直接解析为十进制
    if '°' not in coord_str and '′' not in coord_str and '´' not in coord_str:
        try:
            return float(coord_str)
        except ValueError:
            raise ValueError(f"无法解析坐标: {coord_str}")
    
    # 移除度分符号
    coord_str = coord_str.replace('°', ' ').replace('′', '').replace('´', '')
    parts = coord_str.split()
    
    if len(parts) == 1:
        # 只有度数，没有分
        return float(parts[0])
    elif len(parts) >= 2:
        # 度和分
        degrees = float(parts[0])
        minutes = float(parts[1])
        return degrees + minutes / 60.0
    else:
        raise ValueError(f"无法解析坐标: {coord_str}")

def convert_to_decimal(baseline_data):
    """
    将所有坐标点转换为十进制格式
    """
    coordinates = []
    for name, lat_str, lon_str in baseline_data:
        # 处理北纬
        lat = parse_coordinate(lat_str.replace('北纬', '').replace('N', '').strip())
        # 处理东经
        lon = parse_coordinate(lon_str.replace('东经', '').replace('E', '').strip())
        coordinates.append((name, lat, lon))
    return coordinates

def calculate_polygon_area(coords_latlon):
    """
    计算多边形面积（直线连接方法）
    使用投影坐标系统进行精确计算
    """
    # 创建坐标转换器：从WGS84（经纬度）到UTM Zone 50N（适用于该区域）
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
    
    # 转换坐标到投影坐标系
    projected_coords = []
    for lat, lon in coords_latlon:
        x, y = transformer.transform(lon, lat)
        projected_coords.append((x, y))
    
    # 使用Shapely计算面积
    polygon = Polygon(projected_coords)
    area_m2 = polygon.area
    area_km2 = area_m2 / 1_000_000  # 转换为平方公里
    
    return area_km2, projected_coords

def calculate_distances(coords_latlon):
    """
    计算相邻基线标志点之间的距离，包括闭合回到起点的最后一段
    返回每段距离（单位：公里和海里）
    """
    # 转换到投影坐标系
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
    
    projected_coords = []
    for lat, lon in coords_latlon:
        x, y = transformer.transform(lon, lat)
        projected_coords.append((x, y))
    
    # 计算相邻点之间的距离，包括从最后一点回到第一点
    distances = []
    total_distance_m = 0
    
    for i in range(len(projected_coords)):
        # 当前点
        x1, y1 = projected_coords[i]
        # 下一点（最后一点的下一点是第一点）
        next_idx = (i + 1) % len(projected_coords)
        x2, y2 = projected_coords[next_idx]
        
        distance_m = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        distance_km = distance_m / 1000
        distance_nm = distance_km / 1.852  # 海里换算
        
        distances.append({
            'from_point': i + 1,
            'to_point': next_idx + 1,
            'distance_m': distance_m,
            'distance_km': distance_km,
            'distance_nm': distance_nm
        })
        total_distance_m += distance_m
    
    total_distance_km = total_distance_m / 1000
    total_distance_nm = total_distance_km / 1.852
    
    return distances, total_distance_km, total_distance_nm

def bezier_curve(points, num_points=1000):
    """
    使用贝塞尔曲线插值
    """
    n = len(points) - 1
    t = np.linspace(0, 1, num_points)
    curve = np.zeros((num_points, 2))
    
    for i in range(n + 1):
        bernstein = comb(n, i) * (t ** i) * ((1 - t) ** (n - i))
        curve += np.outer(bernstein, points[i])
    
    return curve[:, 0], curve[:, 1]

def calculate_interpolated_area(coords_latlon, method='cubic_spline', num_points=1000):
    """
    计算使用不同插值方法的曲线围成的面积
    
    参数:
    - method: 插值方法
        - 'linear': 线性插值（直线连接）
        - 'cubic_spline': 三次样条插值
        - 'quadratic_spline': 二次样条插值
        - 'bezier': 贝塞尔曲线
    - num_points: 在曲线上生成的点数
    """
    # 转换到投影坐标系
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
    
    projected_coords = []
    for lat, lon in coords_latlon:
        x, y = transformer.transform(lon, lat)
        projected_coords.append((x, y))
    
    # 提取x和y坐标
    x_coords = [coord[0] for coord in projected_coords]
    y_coords = [coord[1] for coord in projected_coords]
    
    try:
        if method == 'linear':
            # 线性插值（实际上就是直线连接）
            x_new = np.array(x_coords + [x_coords[0]])
            y_new = np.array(y_coords + [y_coords[0]])
            
        elif method == 'cubic_spline':
            # 三次样条插值
            tck, u = splprep([x_coords, y_coords], s=0, k=3, per=True)
            u_new = np.linspace(0, 1, num_points)
            x_new, y_new = splev(u_new, tck)
            
        elif method == 'quadratic_spline':
            # 二次样条插值
            tck, u = splprep([x_coords, y_coords], s=0, k=2, per=True)
            u_new = np.linspace(0, 1, num_points)
            x_new, y_new = splev(u_new, tck)
            
        elif method == 'bezier':
            # 贝塞尔曲线
            points = np.array(projected_coords)
            x_new, y_new = bezier_curve(points, num_points)
            
        else:
            return None, None
        
        # 创建新的多边形并计算面积
        interpolated_coords = list(zip(x_new, y_new))
        interpolated_polygon = Polygon(interpolated_coords)
        area_m2 = interpolated_polygon.area
        area_km2 = area_m2 / 1_000_000
        
        return area_km2, interpolated_coords
    except Exception as e:
        st.error(f"Error calculating {method} interpolation: {e}")
        return None, None

def calculate_smooth_curve_area(coords_latlon, num_points=1000):
    """
    计算光滑曲线围成的面积（三次样条插值）
    为了保持向后兼容
    """
    return calculate_interpolated_area(coords_latlon, method='cubic_spline', num_points=num_points)

def generate_pdf_report(coords_decimal, polygon_area, smooth_area, distances, total_distance_km, 
                        interpolation_results, projected_coords, smooth_coords):
    """
    生成PDF报告，包含可视化结果
    """
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    story = []
    styles = getSampleStyleSheet()
    
    # 标题
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=18,
        textColor=colors.HexColor('#1f77b4'),
        spaceAfter=30,
        alignment=TA_CENTER
    )
    
    story.append(Paragraph("Huangyan Island Territorial Sea Baseline Area Calculation Report", title_style))
    story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles['Normal']))
    story.append(Spacer(1, 20))
    
    # 面积计算结果
    story.append(Paragraph("Area Calculation Results", styles['Heading2']))
    story.append(Spacer(1, 10))
    
    area_data = [
        ['Calculation Method', 'Area (km2)', 'Difference (km2)', 'Percentage (%)'],
        ['Straight-line polygon', f'{polygon_area:.1f}', '-', '-'],
        ['Smooth curve (cubic spline)', f'{smooth_area:.1f}', f'{smooth_area - polygon_area:+.1f}', 
         f'{((smooth_area - polygon_area) / polygon_area * 100):+.2f}']
    ]
    
    area_table = Table(area_data, colWidths=[2.5*inch, 1.5*inch, 1.5*inch, 1.5*inch])
    area_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 10),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))
    story.append(area_table)
    story.append(Spacer(1, 20))
    
    # 周长信息
    story.append(Paragraph("Baseline Perimeter", styles['Heading2']))
    story.append(Paragraph(f"Total perimeter: {total_distance_km:.2f} km ({total_distance_km/1.852:.2f} nm)", styles['Normal']))
    story.append(Spacer(1, 20))
    
    # 可视化对比图
    if projected_coords and smooth_coords:
        story.append(Paragraph("Visualization Comparison", styles['Heading2']))
        story.append(Spacer(1, 10))
        
        try:
            # 生成对比图
            fig = plot_comparison(projected_coords, smooth_coords, coords_decimal)
            img_buffer = BytesIO()
            fig.savefig(img_buffer, format='png', dpi=100, bbox_inches='tight')
            img_buffer.seek(0)
            plt.close(fig)
            
            # 添加图像到PDF
            img = Image(img_buffer, width=7*inch, height=3*inch)
            story.append(img)
            story.append(Spacer(1, 10))
            story.append(Paragraph("Left: Straight-line polygon method. Right: Smooth curve method with baseline points.", 
                                  styles['Normal']))
            story.append(Spacer(1, 20))
        except:
            pass
    
    # 插值方法对比
    if interpolation_results:
        story.append(Paragraph("Interpolation Methods Comparison", styles['Heading2']))
        story.append(Spacer(1, 10))
        
        interp_data = [['Interpolation Method', 'Area (km2)', 'Difference (km2)', 'Percentage (%)']]
        for method_key, result in interpolation_results.items():
            interp_data.append([
                result['name'],
                f"{result['area']:.1f}",
                f"{result['area'] - polygon_area:+.1f}",
                f"{((result['area'] - polygon_area) / polygon_area * 100):+.2f}%"
            ])
        
        interp_table = Table(interp_data, colWidths=[2.5*inch, 1.5*inch, 1.5*inch, 1.5*inch])
        interp_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ]))
        story.append(interp_table)
        story.append(Spacer(1, 20))
    
    # 距离分析
    story.append(Paragraph("Baseline Segment Distances", styles['Heading2']))
    story.append(Spacer(1, 10))
    
    dist_data = [['From Point', 'To Point', 'Distance (km)', 'Distance (nm)']]
    for d in distances:
        dist_data.append([
            f"Point {d['from_point']}",
            f"Point {d['to_point']}",
            f"{d['distance_km']:.3f}",
            f"{d['distance_nm']:.3f}"
        ])
    
    dist_table = Table(dist_data, colWidths=[1.5*inch, 1.5*inch, 1.5*inch, 1.5*inch])
    dist_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTSIZE', (0, 1), (-1, -1), 8)
    ]))
    story.append(dist_table)
    story.append(Spacer(1, 20))
    
    # 基线标志点坐标
    story.append(PageBreak())
    story.append(Paragraph("Baseline Points Coordinates", styles['Heading2']))
    story.append(Spacer(1, 10))
    
    coord_data = [['No.', 'Baseline Point', 'Latitude', 'Longitude']]
    for i, (name, lat, lon) in enumerate(coords_decimal[:-1], 1):
        coord_data.append([str(i), f'Baseline Point {i}', f'{lat:.6f}°', f'{lon:.6f}°'])
    
    coord_table = Table(coord_data, colWidths=[0.5*inch, 2*inch, 2*inch, 2*inch])
    coord_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTSIZE', (0, 1), (-1, -1), 8)
    ]))
    story.append(coord_table)
    
    # 生成PDF
    doc.build(story)
    buffer.seek(0)
    return buffer

def create_interactive_map(coords_decimal):
    """
    创建交互式地图显示基线标志点
    """
    # 计算中心点
    lats = [lat for name, lat, lon in coords_decimal[:-1]]
    lons = [lon for name, lat, lon in coords_decimal[:-1]]
    center_lat = sum(lats) / len(lats)
    center_lon = sum(lons) / len(lons)
    
    # 创建地图
    m = folium.Map(
        location=[center_lat, center_lon],
        zoom_start=12,
        tiles='OpenStreetMap'
    )
    
    # 添加基线标志点
    for i, (name, lat, lon) in enumerate(coords_decimal[:-1]):
        folium.Marker(
            location=[lat, lon],
            popup=f"<b>Baseline Point {i+1}</b><br>Latitude: {lat:.6f}°<br>Longitude: {lon:.6f}°",
            tooltip=f"Point {i+1}",
            icon=folium.Icon(color='red', icon='info-sign')
        ).add_to(m)
    
    # 添加基线多边形
    polygon_coords = [(lat, lon) for name, lat, lon in coords_decimal]
    folium.Polygon(
        locations=polygon_coords,
        color='blue',
        weight=2,
        fill=True,
        fillColor='blue',
        fillOpacity=0.2,
        popup='Huangyan Island Territorial Sea Baseline'
    ).add_to(m)
    
    # 添加图层控制
    folium.LayerControl().add_to(m)
    
    return m

def plot_comparison(projected_coords, smooth_coords, coords_decimal):
    """
    绘制多边形和光滑曲线的对比图
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # 绘制直线连接的多边形
    if projected_coords:
        polygon_array = np.array(projected_coords)
        poly_patch = MplPolygon(polygon_array, fill=True, alpha=0.3, 
                                facecolor='blue', edgecolor='blue', linewidth=2)
        ax1.add_patch(poly_patch)
        ax1.plot(polygon_array[:, 0], polygon_array[:, 1], 'ro-', markersize=8, linewidth=2, label='直线连接')
        ax1.set_title('直线连接方法', fontsize=14, fontweight='bold')
        ax1.set_xlabel('X 坐标 (米)', fontsize=12)
        ax1.set_ylabel('Y 坐标 (米)', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=10)
        ax1.axis('equal')
        
        # 标注点序号
        for i, (name, lat, lon) in enumerate(coords_decimal[:-1]):  # 排除重复的最后一点
            x, y = projected_coords[i]
            ax1.annotate(f'{i+1}', (x, y), xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    # 绘制光滑曲线
    if smooth_coords:
        smooth_array = np.array(smooth_coords)
        poly_patch_smooth = MplPolygon(smooth_array, fill=True, alpha=0.3,
                                       facecolor='green', edgecolor='green', linewidth=2)
        ax2.add_patch(poly_patch_smooth)
        
        # 同时绘制原始点
        if projected_coords:
            polygon_array = np.array(projected_coords)
            ax2.plot(polygon_array[:, 0], polygon_array[:, 1], 'ro', markersize=8, label='基线标志点')
        
        ax2.plot(smooth_array[:, 0], smooth_array[:, 1], 'g-', linewidth=2, label='光滑曲线')
        ax2.set_title('光滑曲线方法', fontsize=14, fontweight='bold')
        ax2.set_xlabel('X 坐标 (米)', fontsize=12)
        ax2.set_ylabel('Y 坐标 (米)', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10)
        ax2.axis('equal')
    
    plt.tight_layout()
    return fig

# 主程序
st.header("📍 领海基线标志点数据")

# 添加自定义数据上传功能
with st.expander("📤 上传自定义基线标志点数据（可选）"):
    st.markdown("""
    您可以上传CSV或Excel文件来计算自定义的基线面积。文件应包含以下列：
    - **标志点名称** (name): 如 "点1", "点2" 等
    - **纬度** (latitude): 可以是度分格式（如 "15°08.1′"）或十进制格式（如 15.135）
    - **经度** (longitude): 可以是度分格式（如 "117°50.9′"）或十进制格式（如 117.848333）
    
    **注意**：最后一个点应与第一个点相同以闭合多边形。
    """)
    
    uploaded_file = st.file_uploader("选择CSV或Excel文件", type=['csv', 'xlsx', 'xls'])
    
    if uploaded_file is not None:
        try:
            # 读取文件
            if uploaded_file.name.endswith('.csv'):
                custom_df = pd.read_csv(uploaded_file)
            else:
                custom_df = pd.read_excel(uploaded_file)
            
            # 检查必需的列
            required_cols = ['name', 'latitude', 'longitude']
            if all(col in custom_df.columns for col in required_cols):
                # 转换为baseline_data格式
                baseline_data = []
                for _, row in custom_df.iterrows():
                    name = str(row['name'])
                    lat_str = str(row['latitude'])
                    lon_str = str(row['longitude'])
                    baseline_data.append((name, lat_str, lon_str))
                
                st.success(f"✅ 成功上传 {len(baseline_data)} 个基线标志点！")
                st.info("使用上传的自定义数据进行计算。")
            else:
                st.error(f"❌ 文件缺少必需的列。需要: {', '.join(required_cols)}")
        except Exception as e:
            st.error(f"❌ 读取文件时出错: {e}")

# 转换坐标
coords_decimal = convert_to_decimal(baseline_data)

# 创建数据表显示
df = pd.DataFrame([
    {
        "序号": i + 1,
        "标志点": name,
        "纬度": f"{lat_str}",
        "经度": f"{lon_str}",
        "纬度(十进制)": f"{coords_decimal[i][1]:.6f}°",
        "经度(十进制)": f"{coords_decimal[i][2]:.6f}°"
    }
    for i, (name, lat_str, lon_str) in enumerate(baseline_data)
])

st.dataframe(df, use_container_width=True, height=600)

# 交互式地图
st.header("🗺️ 交互式地图")
st.markdown("可以缩放、拖动地图查看黄岩岛领海基线的地理位置。点击标记点查看详细坐标信息。")

# 创建并显示交互式地图
interactive_map = create_interactive_map(coords_decimal)
st_folium(interactive_map, width=None, height=500)

# 计算面积
st.header("📐 面积计算结果")

# 提取经纬度坐标（排除重复的最后一点）
coords_latlon = [(lat, lon) for name, lat, lon in coords_decimal[:-1]]

# 计算直线连接的面积
polygon_area, projected_coords = calculate_polygon_area(coords_latlon)

# 计算光滑曲线的面积
smooth_area, smooth_coords = calculate_smooth_curve_area(coords_latlon)

# 显示结果
col1, col2, col3 = st.columns(3)

with col1:
    st.metric(
        label="直线连接方法面积",
        value=f"{polygon_area:.1f} km²",
        help="使用直线连接各基线标志点形成的多边形面积"
    )

with col2:
    if smooth_area is not None:
        st.metric(
            label="光滑曲线方法面积",
            value=f"{smooth_area:.1f} km²",
            help="使用三次样条插值形成的光滑曲线所围面积"
        )

with col3:
    if smooth_area is not None and polygon_area is not None:
        difference = smooth_area - polygon_area
        percentage = (difference / polygon_area) * 100
        st.metric(
            label="面积差异",
            value=f"{abs(difference):.1f} km²",
            delta=f"{percentage:+.2f}%",
            help="光滑曲线面积与直线连接面积的差值"
        )

# 计算基线段距离
distances, total_distance_km, total_distance_nm = calculate_distances(coords_latlon)

# 距离分析
st.header("📏 基线段距离分析")

col1, col2 = st.columns(2)
with col1:
    st.metric(
        label="领海基线总周长",
        value=f"{total_distance_km:.2f} km",
        help="所有基线段的总长度（公里）"
    )
with col2:
    st.metric(
        label="领海基线总周长（海里）",
        value=f"{total_distance_nm:.2f} nm",
        help="所有基线段的总长度（海里，1海里=1.852公里）"
    )

# 显示每段距离的详细信息
with st.expander("📊 查看各基线段详细距离"):
    distance_df = pd.DataFrame([
        {
            "基线段": f"点{d['from_point']} → 点{d['to_point']}",
            "距离(公里)": f"{d['distance_km']:.3f}",
            "距离(海里)": f"{d['distance_nm']:.3f}",
            "距离(米)": f"{d['distance_m']:.1f}"
        }
        for d in distances
    ])
    st.dataframe(distance_df, use_container_width=True)
    
    # 找出最长和最短的基线段
    max_dist = max(distances, key=lambda x: x['distance_km'])
    min_dist = min(distances, key=lambda x: x['distance_km'])
    
    st.info(f"""
    📍 **最长基线段**：点{max_dist['from_point']} → 点{max_dist['to_point']}，长度 {max_dist['distance_km']:.3f} 公里（{max_dist['distance_nm']:.3f} 海里）
    
    📍 **最短基线段**：点{min_dist['from_point']} → 点{min_dist['to_point']}，长度 {min_dist['distance_km']:.3f} 公里（{min_dist['distance_nm']:.3f} 海里）
    """)

# 多种插值方法对比
st.header("🔬 多种插值方法对比分析")

st.markdown("""
对比不同插值方法对面积计算的影响：
- **线性插值**：直线连接各点（即多边形方法）
- **二次样条插值**：使用二次曲线连接
- **三次样条插值**：使用三次曲线连接（默认光滑曲线方法）
- **贝塞尔曲线**：使用贝塞尔曲线插值
""")

# 计算不同插值方法的面积
methods = {
    'linear': 'Linear Interpolation',
    'quadratic_spline': 'Quadratic Spline',
    'cubic_spline': 'Cubic Spline',
    'bezier': 'Bezier Curve'
}

interpolation_results = {}
for method_key, method_name in methods.items():
    area, coords = calculate_interpolated_area(coords_latlon, method=method_key)
    if area is not None:
        interpolation_results[method_key] = {
            'name': method_name,
            'area': area,
            'coords': coords
        }

# 显示对比表格
if interpolation_results:
    comparison_df = pd.DataFrame([
        {
            "插值方法": result['name'],
            "计算面积 (km²)": f"{result['area']:.1f}",
            "与直线方法差异 (km²)": f"{result['area'] - polygon_area:+.1f}",
            "相对差异 (%)": f"{((result['area'] - polygon_area) / polygon_area * 100):+.2f}%"
        }
        for method_key, result in interpolation_results.items()
    ])
    
    st.dataframe(comparison_df, use_container_width=True)
    
    # 可视化对比
    with st.expander("📊 可视化不同插值方法"):
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()
        
        for idx, (method_key, result) in enumerate(interpolation_results.items()):
            ax = axes[idx]
            coords_array = np.array(result['coords'])
            
            # 绘制填充区域
            poly_patch = MplPolygon(coords_array, fill=True, alpha=0.3,
                                   facecolor='blue', edgecolor='blue', linewidth=2)
            ax.add_patch(poly_patch)
            
            # 绘制曲线
            ax.plot(coords_array[:, 0], coords_array[:, 1], 'b-', linewidth=2, label=result['name'])
            
            # 绘制原始标志点
            if projected_coords:
                polygon_array = np.array(projected_coords)
                ax.plot(polygon_array[:, 0], polygon_array[:, 1], 'ro', markersize=8, label='基线标志点')
            
            ax.set_title(f"{result['name']}\n面积: {result['area']:.1f} km²", 
                        fontsize=12, fontweight='bold')
            ax.set_xlabel('X 坐标 (米)', fontsize=10)
            ax.set_ylabel('Y 坐标 (米)', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)
            ax.axis('equal')
        
        plt.tight_layout()
        st.pyplot(fig)
        
    st.info("""
    💡 **分析说明**：
    - 线性插值产生的是多边形，边界为直线段
    - 样条插值方法产生光滑曲线，不同阶数影响曲线的光滑程度
    - 贝塞尔曲线通常会产生更加圆滑的形状，但可能偏离原始点较远
    - 不同方法的面积差异反映了曲线平滑度对结果的影响
    """)

# 详细分析
st.header("📊 详细分析")

if smooth_area is not None and polygon_area is not None:
    difference = smooth_area - polygon_area
    percentage = (difference / polygon_area) * 100
    
    st.markdown(f"""
    ### 计算结果总结：
    
    - **直线连接方法**：将各基线标志点按顺序用直线连接，形成一个多边形。计算得到的面积为 **{polygon_area:.1f} 平方公里**。
    
    - **光滑曲线方法**：使用三次样条插值（Cubic Spline），将基线标志点连接成一条连续光滑的曲线。计算得到的面积为 **{smooth_area:.1f} 平方公里**。
    
    - **面积变化**：光滑曲线方法相比直线连接方法，面积{'增加' if difference > 0 else '减少'}了 **{abs(difference):.1f} 平方公里**，
      相对变化为 **{abs(percentage):.2f}%**。
    
    ### 技术说明：
    
    - 坐标系统：使用 WGS84 (EPSG:4326) 地理坐标系统转换为 UTM Zone 50N (EPSG:32650) 投影坐标系统进行精确面积计算
    - 样条插值：采用三次样条插值（k=3），曲线经过所有基线标志点（s=0），保证闭合（per=True）
    - 面积精度：计算结果精确到 0.1 平方公里
    """)

# 可视化对比
st.header("🗺️ 可视化对比")

if projected_coords and smooth_coords:
    fig = plot_comparison(projected_coords, smooth_coords, coords_decimal)
    st.pyplot(fig)
    
    st.info("💡 左图显示直线连接的多边形，右图显示光滑曲线及基线标志点的位置。")

# 方法说明
with st.expander("📖 计算方法说明"):
    st.markdown("""
    ### 1. 直线连接方法
    
    - 将相邻的基线标志点用直线段连接，形成一个封闭的多边形
    - 使用 Shoelace 公式或 Shapely 几何库计算多边形面积
    - 优点：计算简单，符合实际的领海基线定义
    - 缺点：未考虑可能的曲线特征
    
    ### 2. 光滑曲线方法
    
    - 使用三次样条插值（Cubic Spline）将基线标志点连接成光滑曲线
    - 曲线通过所有标志点，保证连续性和光滑性
    - 在曲线上密集采样（1000个点），形成近似多边形计算面积
    - 优点：考虑了自然地理形态的光滑特征
    - 缺点：可能与实际法律定义的直线基线有偏差
    
    ### 3. 坐标转换
    
    - 经纬度坐标（度分格式）→ 十进制度数
    - WGS84地理坐标系 → UTM投影坐标系（米）
    - UTM Zone 50N 适用于东经 114°-120° 区域，确保面积计算精度
    """)

# PDF导出功能
st.header("📄 导出计算报告")

st.markdown("点击下方按钮生成并下载PDF格式的计算报告，包含所有计算结果和数据表格。")

if st.button("生成PDF报告", type="primary"):
    with st.spinner("正在生成PDF报告..."):
        try:
            pdf_buffer = generate_pdf_report(
                coords_decimal, 
                polygon_area, 
                smooth_area, 
                distances, 
                total_distance_km,
                interpolation_results,
                projected_coords,
                smooth_coords
            )
            
            st.download_button(
                label="📥 下载PDF报告",
                data=pdf_buffer,
                file_name=f"huangyan_island_baseline_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                mime="application/pdf"
            )
            st.success("✅ PDF报告生成成功！点击上方按钮下载。")
        except Exception as e:
            st.error(f"❌ 生成PDF报告时出错: {e}")

st.markdown("---")
st.caption("数据来源：《中华人民共和国政府关于黄岩岛领海基线的声明》（2024年11月10日）")
