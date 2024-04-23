import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import geopandas as gpd
from shapely import Point, MultiPoint, Polygon, LineString, delaunay_triangles, MultiLineString, distance, MultiPolygon, covers, contains_xy, build_area, GeometryCollection
import matplotlib.pyplot as plt
import pandas as pd
import copy
from pyproj import CRS, Transformer
from matplotlib_scalebar.scalebar import ScaleBar
import seaborn as sns
import matplotlib as mpl
import math

# ВСПОМОГАТЕЛЬНЫЕ КЛАССЫ

class Vertice():
    def __init__(self, pnt=(0, 0, 0), i=None, he_list=[], status=None, coast_dist=0) -> None:
        self.pnt = Point(pnt)
        self.i = i
        self.he_list=he_list.copy()
        self.status=status
        self.coast_dist = coast_dist
        self.deep = 0
        self.grad = [None, None]
    def show(self):
        print(self.pnt, self.i, self.he_list, self.status)

class Hedge():
    def __init__(self, i=None, origin=None, face=None, twin=None, next=None, prev=None, edge=None) -> None:
        self.i = i
        self.origin = origin
        self.face = face
        self.twin = twin
        self.next = next
        self.prev = prev
        self.edge = edge
    def show(self):
        print(self.lin, self.i, self.origin, self.face, self.twin, self.next, self.prev)
class Edge():
    def __init__(self, lin=((0, 0), (0, 0)), hedge=None) -> None:
        self.lin = LineString(lin)
        self.hedge=hedge
        self.grad = 0

class Face():
    def __init__(self, pol=(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0)), i=None, p0=None, p1=None, p2=None, outer=None) -> None:
        self.pol = Polygon(pol)
        self.i = i
        self.p0, self.p1, self.p2 = p0, p1, p2
        self.outer = outer
        self.normal = [0, 0, 1]
    def show(self):
        print(self.pol, self.i, self.p0, self.p1, self.p2, self.outer)

# КЛАСС ТРИАНГУЛЯЦИИ

class Triangulation():
    def __init__(self, data, border, crs='EPSG:4326'): # конструктор класса    
        if type(border) == str:
            border = gpd.read_file(border).copy()
        if type(data) == str:
            data = gpd.read_file(data).copy()

        if type(border.geometry[0]) == LineString:
            self.ml = MultiLineString(list(border.geometry))
            mp = MultiPolygon([build_area(GeometryCollection(list(border.geometry)))])
        else: 
            mp = MultiPolygon(list(border.geometry))
            mls = []
            for g in mp.geoms:
                mls.append(g.exterior)
                mls += list(g.interiors)
            self.ml = MultiLineString(mls)
        
        self.crs=crs
        pts_i = []
        self.zs=[]
        self.vertices = []
        self.hedges = []
        self.faces = []
        self.edges = []
        if data.geometry[0] == None:
            for index, i in data.iterrows():
                pts_i.append(Point(float(i.lon), float(i.lat), index))
                self.vertices.append(Vertice((float(i.lon), float(i.lat), float(i.H_sr)), len(self.vertices)))
                self.zs.append(float(i.H_sr))  
        else:
            for index, i in data.iterrows():
                pts_i.append(Point(i.geometry.x, i.geometry.y, index))
                self.vertices.append(Vertice((i.geometry.x, i.geometry.y, float(i.field_3)), len(self.vertices)))
                self.zs.append(float(i.field_3))  
        dT = gpd.GeoSeries(delaunay_triangles(MultiPoint(pts_i)).geoms)

        
        for tr in dT:
            i0 = int(tr.exterior.coords[0][2])
            i1 = int(tr.exterior.coords[1][2])
            i2 = int(tr.exterior.coords[2][2])
            cx = (self.vertices[i0].pnt.x + self.vertices[i1].pnt.x + self.vertices[i2].pnt.x) / 3
            cy = (self.vertices[i0].pnt.y + self.vertices[i1].pnt.y + self.vertices[i2].pnt.y) / 3
            #cross = self.ml.crosses(tr)
            if not contains_xy(mp, cx, cy): #and not cross:
                continue
            self.faces.append(Face((self.vertices[i0].pnt, self.vertices[i1].pnt, self.vertices[i2].pnt, self.vertices[i0].pnt), len(self.faces), i0, i1, i2))
            for p in (i0, i1, i2):
                self.hedges.append(Hedge(len(self.hedges), p, len(self.faces) - 1))
                self.vertices[p].he_list.append(len(self.hedges) - 1)

            e = len(self.hedges)
            self.faces[len(self.faces) - 1].outer = self.hedges[e - 3].i

            eis = [e-1, e-2, e-3]
            for i in range(3):          
                self.hedges[eis[i]].next = eis[i-1]
                self.hedges[eis[i]].prev = eis[i-2]

            for he in range(e - 3, e):
                for i in self.vertices[self.hedges[self.hedges[he].next].origin].he_list:
                    if self.hedges[self.hedges[i].next].origin == self.hedges[he].origin:
                        self.hedges[i].twin = he
                        self.hedges[he].twin = i
                        break

            self.faces[len(self.faces) - 1].outer = e - 3
        for he in self.hedges:
            if he.twin == None:
                he.twin = he.i
        
        he_copy = self.hedges.copy()
        for he in he_copy:
            if he == None:
                continue
            p1 = self.vertices[he.origin].pnt
            p2 = self.vertices[self.hedges[he.next].origin].pnt
            self.edges.append(Edge((Point(p1.x, p1.y), Point(p2.x, p2.y)), he.i))
            self.hedges[he.i].edge = len(self.edges) - 1
            if he.twin != None:
                self.hedges[he.twin].edge = len(self.edges) - 1
                he_copy[he.twin] = None 

    def change_projection(self, new_crs): # функция проецирования
        crs1 = CRS(self.crs)
        crs2 = CRS(new_crs)
        trans = Transformer.from_crs(crs1, crs2)  
        if crs1.is_geographic and crs2.is_projected or crs2.is_geographic and crs1.is_projected:
            for v in self.vertices:
                x, y= trans.transform(v.pnt.y, v.pnt.x, errcheck = True)
                v.pnt = Point(x, y, v.pnt.z)
            self.ml = MultiLineString([LineString([trans.transform(i.xy[1][j], i.xy[0][j], errcheck = True) for j in range(len(i.xy[0]))]) for i in self.ml.geoms])
        else: 
            for v in self.vertices:
                x, y = trans.transform(v.pnt.x, v.pnt.y, errcheck = True)
                v.pnt = Point(x, y, v.pnt.z)
            self.ml = MultiLineString([LineString([trans.transform(i.xy[0][j], i.xy[1][j], errcheck = True) for j in range(len(i.xy[0]))]) for i in self.ml.geoms])
        for f in self.faces:
            f.pol = Polygon((self.vertices[f.p0].pnt,  self.vertices[f.p1].pnt, self.vertices[f.p2].pnt, self.vertices[f.p0].pnt))
        for e in self.edges:
            e.lin = LineString((self.vertices[self.hedges[e.hedge].origin].pnt, self.vertices[self.hedges[self.hedges[e.hedge].next].origin].pnt))
        
        self.crs=new_crs


# БЛОК ФУНКЦИЙ ПОСТРОЕНИЯ


    def calculate_levels(self, levels=[]): # функция построения изолиний и полигонов послойной окраски

        def get_level_pnt(p0, p1, level, linePnts):
            linePnts.append((int(p0.x + (p0.z-level)/(p0.z-p1.z)*(p1.x-p0.x)), int(p0.y + (p0.z-level)/(p0.z-p1.z)*(p1.y-p0.y))))

        def get_level_line(tr, level, trs_viewed):
            linePnts = []
            #print("New_line!")
            p0 = self.vertices[self.hedges[tr.outer].origin].pnt
            p1 = self.vertices[self.hedges[self.hedges[tr.outer].next].origin].pnt
            p2 = self.vertices[self.hedges[self.hedges[tr.outer].prev].origin].pnt
            if (p0.z - level) * (p1.z - level) < 0 or (p0.z - level) * (p2.z - level) < 0 or (p2.z - level) * (p1.z - level) < 0:
                trs_viewed[tr.i] = 1
                ps = (p0, p1, p2)
                hes = {0: self.hedges[self.hedges[tr.outer].prev].twin, 1:self.hedges[tr.outer].twin, 2:self.hedges[self.hedges[tr.outer].next].twin}
                he_i = None
                for i in hes:
                    if (ps[i].z - level) * (ps[i-2].z - level) >= 0:
                        get_level_pnt(ps[i-2], ps[i-1], level, linePnts)
                        get_level_pnt(ps[i], ps[i-1], level, linePnts)
                        he_i = hes[i]
                        break
                if he_i == None:
                    print("Первый иф дал сбой!!!")
                    return "Error"
                #print(tr.i, he_i)

                afterboarder = False
                start = True
                n = 0
                while self.hedges[he_i].face != tr.i and n < len(self.faces) or start:
                    n+=1
                    start = False
                    if he_i == self.hedges[he_i].twin and not afterboarder:
                        #print('aaa')
                        if self.vertices[self.hedges[he_i].origin].pnt.z > level:
                            p_next = self.vertices[self.hedges[he_i].origin]
                        elif self.vertices[self.hedges[self.hedges[he_i].next].origin].pnt.z > level:
                            p_next = self.vertices[self.hedges[self.hedges[he_i].next].origin]
                        else:
                            print("У последнего треугольника горизонталь не на ребре!!!")
                            return "Error"
                            
                        he_i_next = he_i
                        while p_next.pnt.z > level:
                            p = p_next
                            he_i = he_i_next
                            linePnts.append((p.pnt.x, p.pnt.y))
                            #print(f'making a: {self.hedges[he_i].face}, {he_i}, {p.i}')
                            
                            for he in p.he_list:
                                if he == self.hedges[he].twin and he != he_i:
                                    he_i_next = he
                                    p_next = self.vertices[self.hedges[self.hedges[he_i_next].next].origin]
                                    break
                                elif self.hedges[he].prev == self.hedges[self.hedges[he].prev].twin and self.hedges[he].prev != he_i:
                                    he_i_next = self.hedges[he].prev
                                    p_next = self.vertices[self.hedges[he_i_next].origin]
                                    break
                                else:
                                    continue
                        he_i = he_i_next
                        get_level_pnt(p_next.pnt, p.pnt, level, linePnts)
                        #print(f'end a: {self.hedges[he_i].face}, {he_i}')
                        afterboarder = True
                    else:
                        #print(self.hedges[he_i].face, he_i)
                        trs_viewed[self.hedges[he_i].face] = 1
                        p0 = self.vertices[self.hedges[he_i].origin].pnt
                        p1 = self.vertices[self.hedges[self.hedges[he_i].next].origin].pnt
                        p2 = self.vertices[self.hedges[self.hedges[he_i].prev].origin].pnt
                        if (p0.z - level) * (p2.z - level) >= 0: 
                            get_level_pnt(p1, p2, level, linePnts)
                            he_i = self.hedges[self.hedges[he_i].next].twin
                        elif (p1.z - level) * (p2.z - level) >= 0:
                            get_level_pnt(p0, p2, level, linePnts)
                            he_i = self.hedges[self.hedges[he_i].prev].twin
                        else:
                            print("У следующего треугольника нет горизонтали!!!")
                            ax = gpd.GeoSeries([v.pnt for v in self.vertices]).plot(figsize=[200, 100])
                            gpd.GeoSeries.plot(gpd.GeoSeries(LineString(linePnts)), ax=ax)
                            gpd.GeoSeries.plot(gpd.GeoSeries([i.lin for i in self.edges]), ax=ax)
                            plt.draw()
                            return "Error"
                                                  
                        afterboarder = False

                if n == len(self.faces):
                    print("Бесконечный цикл!!!")
                    ax = gpd.GeoSeries([v.pnt for v in self.vertices]).plot(figsize=[200, 100])
                    gpd.GeoSeries.plot(gpd.GeoSeries(LineString(linePnts)), ax=ax)
                    gpd.GeoSeries.plot(gpd.GeoSeries(Point(linePnts[-1])), ax=ax, color='r')
                    plt.draw()
                    return "Error"
                
                return LineString(linePnts)
            else:
                return None
            
        self.contours = dict()
        self.areas = dict()
        self.levels = levels
        for v in self.vertices:
            if v.pnt.z in levels:
                v.pnt = Point(v.pnt.x, v.pnt.y, v.pnt.z+0.001)
        for level in levels:
            trs_viewed = [0] * len(self.faces)
            self.contours[level] = []
            for tr in self.faces:
                if trs_viewed[tr.i] == 1:
                    continue
                line = get_level_line(tr, level, trs_viewed)
                if line == "Error":
                    print("\n\n\nError\n\n\n")
                    break
                if line != None:
                    self.contours[level].append(line)
                    if not line.is_closed:
                        print(line)
              
            print(GeometryCollection(self.contours[level]))
            self.areas[level] = build_area(GeometryCollection(self.contours[level]))
            for i in self.ml.geoms:
                if covers(self.areas[level], i):
                    self.areas[level] = build_area(GeometryCollection([self.areas[level], i]))
            self.areas[level] = self.areas[level] if self.contours[level] != [] else build_area(GeometryCollection(self.ml))
            self.contours[level] = MultiLineString(self.contours[level])
        
    def calculate_distance(self, coast): # функция расчёта значений расстояния до берега
        t = copy.deepcopy(self)
        for v in self.vertices:
            v.coast_dist = distance(coast, Point([v.pnt.x, v.pnt.y]))
            t.vertices[v.i].pnt = Point(v.pnt.x, v.pnt.y, v.coast_dist/1000000) #в тысячах км
        t.calculate_grad()
        for v in self.vertices:
            v.grad_dist = t.vertices[v.i].grad
        
    def calculate_grad(self): # функция расчёта градиентов поля характеристик
        for tr in self.faces:
            p0 = self.vertices[tr.p0].pnt
            p1 = self.vertices[tr.p1].pnt
            p2 = self.vertices[tr.p2].pnt
            longx = (p1.y - p0.y) * (p2.z - p0.z) - (p2.y - p0.y) * (p1.z - p0.z)
            longy = (p1.z - p0.z) * (p2.x - p0.x) - (p2.z - p0.z) * (p1.x - p0.x)
            longz = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)
            leng = math.sqrt(longx**2 + longy**2 + longz**2)
            tr.normal = [longx/leng, longy/leng, longz/leng]
        for e in self.edges:
            e.grad = (self.vertices[self.hedges[e.hedge].origin].pnt.z - self.vertices[self.hedges[self.hedges[e.hedge].next].origin].pnt.z) / e.lin.length
        for v in self.vertices:
            grad = [0, 0, 0]
            for he in v.he_list:
                for i in range(3):
                    grad[i] += self.faces[self.hedges[he].face].normal[i] * self.faces[self.hedges[he].face].pol.area / (
                        self.edges[self.hedges[he].edge].lin.length**2 * self.edges[self.hedges[self.hedges[he].prev].edge].lin.length**2)
            leng = math.sqrt(grad[0]**2 + grad[1]**2 + grad[2]**2)
            leng = leng if leng != 0 else 1
            v.grad = [grad[0]/leng, grad[1]/leng]


# БЛОК ФУНКЦИЙ ГЕНЕРАЛИЗАЦИИ


    def push_value_by_min(self): # функция утрирования значений триангуляции методом ближайшего к берегу соседа
        self.calculate_distance()
        new_t = copy.deepcopy(self)
        for i in range(len(new_t.vertices)):
            z = self.vertices[i].pnt.z
            min = self.vertices[i].coast_dist
            for he in self.vertices[i].he_list:
                if self.vertices[self.hedges[self.hedges[he].next].origin].coast_dist > min:
                    continue
                min = self.vertices[self.hedges[self.hedges[he].next].origin].coast_dist
                z = self.vertices[self.hedges[self.hedges[he].next].origin].pnt.z

            new_t.vertices[i].pnt = Point((self.vertices[i].pnt.x, self.vertices[i].pnt.y, z))
        new_t.zs = [i.pnt.z for i in new_t.vertices]
        return new_t
    
    def push_value_by_dist(self, coast, q=0.5): # функция утрирования значений триангуляции методом подсчёта расстояния до берега
        self.calculate_distance(coast)
        new_t = copy.deepcopy(self)
        max_dist = max([v.coast_dist for v in self.vertices])
        lengths = sorted([e.lin.length for e in self.edges])
        maxl = max(lengths)
        minl = min(lengths)
        norml = maxl - minl

        for new_v in new_t.vertices:
            p0 = self.vertices[new_v.i].pnt
            g0 = (self.vertices[new_v.i].grad_dist[0], self.vertices[new_v.i].grad_dist[1])
            suma = 0
            sumdzi = 0
            kd = 1 - (new_v.coast_dist / max_dist)**q
            for he in self.vertices[new_v.i].he_list:
                pi = self.vertices[self.hedges[self.hedges[he].next].origin].pnt
                vi = (pi.x - p0.x, pi.y - p0.y)
                dot = vi[0]*g0[0] + vi[1]*g0[1]
                dzi = pi.z - p0.z
                if dot <= 0 or dzi > 0:
                     continue
                ei = self.edges[self.hedges[he].edge]
                
                kl = 1 - (ei.lin.length - minl) / norml
                ka = dot / (math.sqrt(g0[0]**2+g0[1]**2) * ei.lin.length)
                
                sumdzi += dzi * kl * ka
                suma += ka
            suma = suma if suma != 0 else 1
            dz = sumdzi * kd / suma
            new_z = p0.z + dz

            new_v.pnt = Point((new_v.pnt.x, new_v.pnt.y, new_z))
        new_t.zs = [i.pnt.z for i in new_t.vertices]
        return new_t

    def push_value_by_grad(self, s): # функция утрирования значений триангуляции методом градиента
        self.calculate_grad()
        new_t = copy.deepcopy(self)
        grads = sorted([abs(e.grad) for e in self.edges])
        normg = grads[len(grads) - len(grads)//100]
        
        lengths = sorted([e.lin.length for e in self.edges])
        maxl = max(lengths)
        minl = min(lengths)
        norml = maxl - minl

        for new_v in new_t.vertices:
            p0 = self.vertices[new_v.i].pnt
            g0 = (self.vertices[new_v.i].grad[0], self.vertices[new_v.i].grad[1])
            suma = 0
            sumdzi = 0
            for he in self.vertices[new_v.i].he_list:
                pi = self.vertices[self.hedges[self.hedges[he].next].origin].pnt
                vi = (pi.x - p0.x, pi.y - p0.y)
                dot = vi[0]*g0[0] + vi[1]*g0[1]
                dzi = pi.z - p0.z
                if dot <= 0 or dzi > 0:
                     continue
                ei = self.edges[self.hedges[he].edge]
                
                kl = 1 - (ei.lin.length - minl) / norml
                kg = (abs(ei.grad) / normg)**s if abs(ei.grad) / normg < 1 else 1
                ka = dot / (math.sqrt(g0[0]**2+g0[1]**2) * ei.lin.length)
                
                sumdzi += dzi * 1 * kg * ka
                suma += ka
            suma = suma if suma != 0 else 1
            dz = sumdzi / suma
            new_z = p0.z + dz

            new_v.pnt = Point((new_v.pnt.x, new_v.pnt.y, new_z))
        new_t.zs = [i.pnt.z for i in new_t.vertices]
        return new_t

    def smoother(self): # функция сглаживания триангуляции методом среднего арифметического
        new_t = copy.deepcopy(self)
        for v in new_t.vertices:
            s = self.vertices[v.i].pnt.z
            n = 1
            for he in self.vertices[v.i].he_list:
                s += self.vertices[self.hedges[self.hedges[he].next].origin].pnt.z
                if self.hedges[self.hedges[he].prev].i == self.hedges[self.hedges[he].prev].twin:
                    s += self.vertices[self.hedges[self.hedges[he].prev].origin].pnt.z
                    n += 1    
                n += 1
            v.pnt = Point((self.vertices[v.i].pnt.x, self.vertices[v.i].pnt.y, s / n))
        new_t.zs = [i.pnt.z for i in new_t.vertices]
        return new_t

    def smart_smoother(self, depth=1): # функция сглаживания триангуляции методом ОВР
        def next_neighbor(vi, found_set, depth, max_depth, c):
            dist = math.sqrt((c[0] - self.vertices[vi].pnt.x)**2 + (c[1] - self.vertices[vi].pnt.y)**2)
            if depth == max_depth:
                 k = (1 / dist)
                 return self.vertices[vi].pnt.z * k, k, set([vi])
            v_set = set()
            for he in self.vertices[vi].he_list:
                 v_next = self.hedges[self.hedges[he].next].origin
                 if self.hedges[self.hedges[he].prev].i == self.hedges[self.hedges[he].prev].twin and not (self.hedges[self.hedges[he].prev].origin in found_set):
                    v_set.add(self.hedges[self.hedges[he].prev].origin)
                 if not (v_next in found_set):
                    v_set.add(v_next) 
            v_sumk = (1 / dist) if dist > 0 else 0
            v_sumz = self.vertices[vi].pnt.z * v_sumk
            v_sumset = v_set
            for v in v_set:
                new_z, new_k, new_vs = next_neighbor(v, found_set | v_sumset, depth + 1, max_depth, c)
                v_sumz += new_z
                v_sumk += new_k
                v_sumset = new_vs | v_sumset
            return v_sumz, v_sumk, v_sumset
                             
        new_t = copy.deepcopy(self)
        for v in new_t.vertices:
            s, k, neighbours_set = next_neighbor(v.i, set([v.i]), 0, max_depth=depth, c=[v.pnt.x, v.pnt.y])
            s *= 1/k
            v.pnt = Point((v.pnt.x, v.pnt.y, (s + self.vertices[v.i].pnt.z)/ 2))
        new_t.zs = [i.pnt.z for i in new_t.vertices]
        return new_t


# БЛОК ФУНКЦИЙ ВИЗУАЛИЗАЦИИ


    def show(self, ax):
        gpd.GeoSeries.plot(gpd.GeoSeries([i.lin for i in self.edges]), ax=ax)
        gpd.GeoSeries.plot(gpd.GeoSeries([p.pnt for p in self.vertices]), ax=ax)

    def show_coast(self, ax):
        gpd.GeoSeries.plot(gpd.GeoSeries(self.ml), ax=ax)
        gpd.GeoSeries.plot(gpd.GeoSeries([p.pnt for p in self.vertices]), ax=ax)

    def draw_levels(self, ax, ground, sea, cmap='BuGn'):
        #gpd.GeoSeries.plot(gpd.GeoSeries([p.pnt for p in self.vertices]), ax = ax)
        gpd.GeoSeries.plot(gpd.GeoSeries([self.areas[i] for i in self.levels]), cmap=cmap, ax=ax, zorder=1)
        gpd.GeoSeries.plot(gpd.GeoSeries([self.contours[i] for i in self.contours]), color = "k", ax=ax, linewidth = 1, zorder =2)
        gpd.GeoSeries.plot(sea, ax=ax, color='#f1f1f5', edgecolor='C0', linewidth = 1, zorder = 4)
        gpd.GeoSeries.plot(ground, ax=ax, color='#afafaf', edgecolor='C0', linewidth = 1, zorder = 5)
        gpd.GeoSeries.plot(gpd.GeoSeries([self.ml]), ax=ax, color='C0', linewidth = 1, zorder = 3)

    def show_grad(self, ax):
        t1 = copy.deepcopy(self)
        for p in t1.vertices:
            p.pnt = Point(p.pnt.x, p.pnt.y, math.sqrt(p.grad[0]**2+p.grad[1]**2))
        gradlist = [math.sqrt(p.grad[0]**2+p.grad[1]**2) for p in t1.vertices]
        self.levels = [0, (min(gradlist)+max(gradlist)) / 2]
        t1.calculate_levels(levels=self.levels)
        print(self.levels)
        t1.draw_levels(ax)


# БЛОК ФУНКЦИЙ СОХРАНЕНИЯ МОДЕЛИ


    def save_to_shp(self, path): # функция сохранения триангуляции в формат shp
        d = {'par': [i.pnt.z for i in self.vertices], 'geometry': [i.pnt for i in self.vertices]}
        gpd.GeoDataFrame(d, crs=self.crs).to_file(path)
        
    def save_to_csv(self, path): # функция сохранения триангуляции в формат csv
        xes1 = pd.Series([i.pnt.x for i in self.vertices])
        yes1 = pd.Series([i.pnt.y for i in self.vertices])
        zes1 = pd.Series([i.pnt.z for i in self.vertices])
        pd.DataFrame(data={'x':xes1, 'y':yes1, 'z':zes1}).to_csv(path)

# КОНЕЦ КЛАССА ТРИАНГУЛЯЦИИ

def make_scaleImg_series(t, xc, yc, ground, sea, k=1, mode='izo', cmap='BuGn', name="Средняя высота волны, м"): #функция визуализации списка из 4 триангуляций в разных масштабах
    fig, ax = plt.subplots(2, 2, figsize=[10, 10])
    move = [90000 * k, 200000 * k, 600000 * k, 1100000 * k]

    for i in range(4):
        ax[i//2, i%2].set_xlim(xc - move[i], xc + move[i])
        ax[i//2, i%2].set_ylim(yc - move[i], yc + move[i])
        ax[i//2, i%2].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax[i//2, i%2].add_artist(ScaleBar(dx=1, location='lower right'))
        if mode == 'izo':
            t[i].draw_levels(ax[i//2, i%2], ground=ground, sea=sea, cmap=cmap)
        else:
            t[i].show(ax[i//2, i%2])
        #plt.subplots_adjust(wspace=0.01, hspace=0.01)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots(figsize=(3, 0.3))
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.BoundaryNorm(t[0].levels + [max(t[0].zs)], mpl.cm.Reds.N), cmap=cmap),
                  label=name, drawedges=True, format="{x:.1f}", orientation='horizontal', cax=ax, shrink=0.5)
    


def make_closeImg_series(t, xc, yc, ground, sea, ext=90000, cmap='BuGn', name="Средняя высота волны, м"): #функция визуализации списка из 4 триангуляций в одном масштабе
    fig, ax = plt.subplots(2, 2, figsize=[10, 10])
    move = [ext] * 4
    titles = []
  
    for i in range(4):
        ax[i//2, i%2].set_xlim(xc - move[i], xc + move[i])
        ax[i//2, i%2].set_ylim(yc - move[i], yc + move[i])
        ax[i//2, i%2].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax[i//2, i%2].add_artist(ScaleBar(dx=1, location='lower right'))
        t[i].draw_levels(ax[i//2, i%2], ground=ground, sea=sea, cmap=cmap)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 0.3))
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.BoundaryNorm(t[0].levels + [max(t[0].zs)], mpl.cm.Reds.N), cmap=cmap),
                  label=name, drawedges=True, format="{x:.5f}", orientation='horizontal', cax=ax, shrink=0.5)
    


def draw_ecdf(data, ax, c = 'b'): # функция построения эмпирической функции распределения
    data.sort()
    y = []
    x = np.linspace(min(data), max(data), 1000)
    j_last = 0
    for i in x:
        for j in range(j_last, len(data)):
            if data[j] >= i:
                y.append(j/len(data))
                j_last = j
                break
            else:
                pass
    print(y)
    ax.plot(x, y, c=c)
def draw_ecdf_zs(ts): # функции построения ЭФР характеристик для списка триангуляций на одном графике
    fig, ax = plt.subplots(figsize=[7, 7])
    colors = sns.color_palette("Reds", len(ts))
    for i in range(len(ts)):
        grads = ts[i].zs
        draw_ecdf(grads, ax, c = colors[i])
def draw_ecdf_grad(ts): # функции построения ЭФР градиента характеристик для списка триангуляций на одном графике
    fig, ax = plt.subplots(figsize=[7, 7])
    colors = sns.color_palette("Reds", len(ts))
    for i in range(len(ts)):
        grads = [abs(e.grad) for e in ts[i].edges]
        draw_ecdf(grads, ax, c = colors[i])