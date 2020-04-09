#include <mbgl/style/expression/distance.hpp>

#include <mapbox/cheap_ruler.hpp>
#include <mapbox/geojson.hpp>
#include <mapbox/geometry.hpp>

#include <mbgl/style/conversion/json.hpp>
#include <mbgl/tile/geometry_tile_data.hpp>
#include <mbgl/util/geometry_within.hpp>

#include <mbgl/util/logging.hpp>
#include <mbgl/util/string.hpp>

#include <rapidjson/document.h>

#include <limits>
#include <queue>
#include <tuple>

namespace mbgl {
namespace {
//
// using namespace mapbox::geometry;
//
// void updateBBox(std::array<double, 4>& bbox, point<double> p) {
//       bbox[0] = std::min(p.x, bbox[0]);
//       bbox[1] = std::min(p.y, bbox[1]);
//       bbox[2] = std::max(p.x, bbox[2]);
//       bbox[3] = std::max(p.y, bbox[3]);
//}
//
// std::array<double, 4> getBBox( const multi_point<double>& points){
//    std::array<double, 4> bbox{std::numeric_limits<double>::infinity(),
//                                                        std::numeric_limits<double>::infinity(),
//                                                        -std::numeric_limits<double>::infinity(),
//                                                       -std::numeric_limits<double>::infinity()};
//    for(const auto& p: points){
//        updateBBox(bbox, p);
//    }
//    return bbox;
//}
//// Change to use unit
//// bbox[minX, minY, maxX, maxY]
// double calculateBBoxDistance(std::array<double, 4>& bbox1, std::array<double, 4>& bbox2,
// mapbox::cheap_ruler::CheapRuler& ruler) {
//
//    // bbox1 in left side
//    if(bbox1[2] < bbox2[0]) {
//        if(bbox1[1] > bbox2[3]) {
//            return ruler.distance(point<double>{bbox1[2], bbox1[1]}, point<double>{bbox2[0], bbox2[3]});
//        }
//        if(bbox1[3] < bbox2[1]) {
//            return ruler.distance(point<double>{bbox1[2], bbox1[3]}, point<double>{bbox2[0], bbox2[1]});
//        }
//        auto p1 = point<double>{bbox1[2], bbox1[1]};
//        auto p2 = point<double>{bbox1[2], bbox1[3]};
//        auto q1 = point<double>{bbox2[0], bbox2[3]};
//        auto q2 = point<double>{bbox2[0], bbox2[1]};
//        auto min = std::min(ruler.distance(p1, q1), ruler.distance(p1, q2));
//        min = std::min(min, ruler.distance(p2, q1));
//        min = std::min(min, ruler.distance(p2, q2));
//        return min;
//    }
//
//    // bbox1 in right side
//    if(bbox1[0] > bbox2[2]) {
//        if(bbox1[1] > bbox2[3]) {
//           return ruler.distance(point<double>{bbox1[0], bbox1[1]}, point<double>{bbox2[2], bbox2[3]});
//       }
//       if(bbox1[3] < bbox2[1]) {
//           return ruler.distance(point<double>{bbox1[0], bbox1[3]}, point<double>{bbox2[2], bbox2[1]});
//       }
//       auto p1 = point<double>{bbox1[0], bbox1[1]};
//       auto p2 = point<double>{bbox1[0], bbox1[3]};
//       auto q1 = point<double>{bbox2[2], bbox2[3]};
//       auto q2 = point<double>{bbox2[2], bbox2[1]};
//       auto min = std::min(ruler.distance(p1, q1), ruler.distance(p1, q2));
//       min = std::min(min, ruler.distance(p2, q1));
//       min = std::min(min, ruler.distance(p2, q2));
//       return min;
//    }
//    double min =  std::numeric_limits<double>::infinity();
//    for(std::size_t i = 0; i <= ) {
//        for(const auto& q: bbox2) {
//            min = std::min(min, ruler.distance(p, q));
//        }
//    }
//
//    return min;
//}

// https://stackoverflow.com/questions/45861488/distance-between-two-polylines
// Divide and conqure, the time complexity is O(n*lgn), faster than O(n*n)
// however it requires extra space, trade off
// double pointsToPointsDistance(const mapbox::geometry::multi_point<double>& points1,
//                                const mapbox::geometry::multi_point<double>& points2,
//                                mapbox::cheap_ruler::CheapRuler& ruler) {
//    using DistPair = std::tuple<double, multi_point<double>, multi_point<double>>;
//    using DistQueue = std::queue<DistPair>;
//
//    auto miniDist = ruler.distance(points1[0], points2[0]);
//
//    DistQueue distQueue;
//    distQueue.push(std::make_tuple(0, points1, points2));
//
//    while(!distQueue.empty()) {
//        const auto& distPair = distQueue.front();
//        distQueue.pop();
//        if (std::get<0>(distPair) >= miniDist) break;
//        auto& pSetA = std::get<1>(distPair);
//        auto& pSetB = std::get<2>(distPair);
//        static const std::size_t MinClusterSize = 5;
//        if (pSetA.size() <= MinClusterSize ||  pSetB.size() <= MinClusterSize) {
//            for (const auto& p: pSetA) {
//                auto temp = pointToPointsDistance(p, pSetB, ruler);
//                if (temp == 0.) {
//                    miniDist = 0.;
//                    break;
//                }
//                miniDist = std::min(miniDist, temp);
//            }
//        } else {
//            auto size1 = pSetA.size() / 2;
//            auto pSetA1 = multi_point<double>(pSetA.begin(),  pSetA.begin() + size1);
//            auto pSetA2 = multi_point<double>( pSetA.begin() + size1,  pSetA.end());
////            pSetA.clear();
//
//            auto size2 = pSetB.size() / 2;
//            auto pSetB1 = multi_point<double>(pSetB.begin(),  pSetB.begin() + size2);
//            auto pSetB2 = multi_point<double>( pSetB.begin() + size2,  pSetB.end());
////            pSetB.clear();
//
//            const auto updateQueue = [&distQueue, &miniDist](const multi_point<double>& set1, const
//            multi_point<double>& set2) {
//
//                auto tempDist = calculateBBoxDistance(getBBox(set1), getBBox(set2));
//                if (tempDist <= miniDist) insertQueue(distQueue, std::make_tuple(tempDist, set1, set2));
//            };
//
//            updateQueue(pSetA1, pSetB1);
//            updateQueue(pSetA1, pSetB2);
//            updateQueue(pSetA2, pSetB1);
//            updateQueue(pSetA2, pSetB2);
//
//        }
//
//    }
//    return miniDist;
//}

double pointToLineDistance(const mapbox::geometry::point<double>& point,
                           const mapbox::geometry::line_string<double>& line,
                           mapbox::cheap_ruler::CheapRuler& ruler) {
    const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
    return ruler.distance(point, nearestPoint);
}

double pointToLinesDistance(const mapbox::geometry::point<double>& point,
                            const mapbox::geometry::multi_line_string<double>& lines,
                            mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& line : lines) {
        const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
        auto tempDist = ruler.distance(point, nearestPoint);
        if (tempDist == 0.) return tempDist;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double pointToPointsDistance(const mapbox::geometry::point<double>& point,
                             const mapbox::geometry::multi_point<double>& points,
                             mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < points.size(); ++i) {
        auto tempDist = ruler.distance(point, points[i]);
        if (tempDist == 0.) return tempDist;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double pointToPolygonDistance(const mapbox::geometry::point<double>& point,
                              const mapbox::geometry::polygon<double>& polygon,
                              mapbox::cheap_ruler::CheapRuler& ruler) {
    // check if point is within polygon or not
    if (GeometryUtil<double>().pointWithinPolygon(point, polygon, true /*checkBoundary*/)) return 0.;

    double dist = std::numeric_limits<double>::infinity();
    for (const auto& line : polygon) {
        const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
        auto tempDist = ruler.distance(point, nearestPoint);
        if (tempDist == 0.) return tempDist;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double pointToPolygonsDistance(const mapbox::geometry::point<double>& point,
                               const mapbox::geometry::multi_polygon<double>& polygons,
                               mapbox::cheap_ruler::CheapRuler& ruler) {
    // check if point is within polygon or not
    for (const auto& polygon : polygons) {
        if (GeometryUtil<double>().pointWithinPolygon(point, polygon, true /*checkBoundary*/)) return 0.;
    }

    double dist = std::numeric_limits<double>::infinity();
    for (const auto& polygon : polygons) {
        for (const auto& line : polygon) {
            const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
            auto tempDist = ruler.distance(point, nearestPoint);
            if (tempDist == 0.) return tempDist;
            dist = std::min(dist, tempDist);
        }
    }
    return dist;
}

double lineToLineDistance(const mapbox::geometry::line_string<double>& line1,
                          const mapbox::geometry::line_string<double>& line2,
                          mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < line1.size() - 1; ++i) {
        const auto& p1 = line1[i];
        const auto& p2 = line1[i + 1];
        for (std::size_t j = 0; j < line2.size() - 1; ++j) {
            const auto& q1 = line2[j];
            const auto& q2 = line2[j + 1];

            if (GeometryUtil<double>().segmentIntersectSegment(p1, p2, q1, q2)) return 0.;
            dist = std::min(dist, pointToLineDistance(p1, mapbox::geometry::line_string<double>{q1, q2}, ruler));
            dist = std::min(dist, pointToLineDistance(p2, mapbox::geometry::line_string<double>{q1, q2}, ruler));
            dist = std::min(dist, pointToLineDistance(q1, mapbox::geometry::line_string<double>{p1, p2}, ruler));
            dist = std::min(dist, pointToLineDistance(q1, mapbox::geometry::line_string<double>{p1, p2}, ruler));
        }
    }
    return dist;
}

double lineToLinesDistance(const mapbox::geometry::line_string<double>& line,
                           const mapbox::geometry::multi_line_string<double>& lines,
                           mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& l : lines) {
        auto tempDist = lineToLineDistance(line, l, ruler);
        if (tempDist == 0.) return 0.;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double lineToPolygonDistance(const mapbox::geometry::line_string<double>& line,
                           const mapbox::geometry::polygon<double>& polygon,
                           mapbox::cheap_ruler::CheapRuler& ruler) {
    double dist = std::numeric_limits<double>::infinity();
    for (const auto& l : polygon) {
        auto tempDist = lineToLineDistance(line, l, ruler);
        if (tempDist == 0.) return 0.;
        dist = std::min(dist, tempDist);
    }
    return dist;
}

double pointToGeometryDistance(const mapbox::geometry::point<double>& point,
                               const Feature::geometry_type& geoSet,
                               mapbox::cheap_ruler::CheapRuler::Unit unit) {
    mapbox::cheap_ruler::CheapRuler ruler(point.y, unit);
    return geoSet.match([&point, &ruler](const mapbox::geometry::point<double>& p) { return ruler.distance(point, p); },
                        [&point, &ruler](const mapbox::geometry::multi_point<double>& points) {
                            return pointToPointsDistance(point, points, ruler);
                        },
                        [&point, &ruler](const mapbox::geometry::line_string<double>& line) {
                            return pointToLineDistance(point, line, ruler);
                        },
                        [&point, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
                            return pointToLinesDistance(point, lines, ruler);
                        },
                        [&point, &ruler](const mapbox::geometry::polygon<double>& polygon) {
                            return pointToPolygonDistance(point, polygon, ruler);
                        },
                        [&point, &ruler](const mapbox::geometry::multi_polygon<double>& polygons) {
                            return pointToPolygonsDistance(point, polygons, ruler);
                        },
                        [](const auto&) -> double { return std::numeric_limits<double>::quiet_NaN(); });
}

double lineToGeometryDistance(const mapbox::geometry::line_string<double>& line,
                              const Feature::geometry_type& geoSet,
                              mapbox::cheap_ruler::CheapRuler::Unit unit) {
    assert(!line.empty());
    mapbox::cheap_ruler::CheapRuler ruler(line.front().y, unit);
    return geoSet.match(
        [&line, &ruler](const mapbox::geometry::point<double>& p) { return pointToLineDistance(p, line, ruler); },
        [&line, &ruler](const mapbox::geometry::multi_point<double>& points) {
            double dist = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < points.size(); ++i) {
                auto tempDist = pointToLineDistance(points[i], line, ruler);
                if (tempDist == 0.) return 0.;
                dist = std::min(dist, tempDist);
            }
            return dist;
        },
        [&line, &ruler](const mapbox::geometry::line_string<double>& line1) {
            return lineToLineDistance(line, line1, ruler);
        },
        [&line, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
            return lineToLinesDistance(line, lines, ruler);
        },
        [&line, &ruler](const mapbox::geometry::polygon<double>& polygon) {
            return lineToPolygonDistance(line, polygon, ruler);
        },
        [&line, &ruler](const mapbox::geometry::multi_polygon<double>& polygons) {
            double dist = std::numeric_limits<double>::infinity();
            for (const auto& polygon : polygons) {
                auto tempDist = lineToPolygonDistance(line, polygon, ruler);
                if (tempDist == 0.) return 0.;
                dist = std::min(dist, tempDist);
            }
            return dist;
        },
        [](const auto&) -> double { return std::numeric_limits<double>::quiet_NaN(); });
}

double polygonToGeometryDistance(const mapbox::geometry::polygon<double>& polygon,
                                 const Feature::geometry_type& geoSet,
                                 mapbox::cheap_ruler::CheapRuler::Unit unit) {
    assert(!polygon.empty());
    mapbox::cheap_ruler::CheapRuler ruler(polygon.front().front().y, unit);
    return geoSet.match(
        [&polygon, &ruler](const mapbox::geometry::point<double>& p) {
            return pointToPolygonDistance(p, polygon, ruler);
        },
        [&polygon, &ruler](const mapbox::geometry::multi_point<double>& points) {
            double dist = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < points.size(); ++i) {
                auto tempDist = pointToPolygonDistance(points[i], polygon, ruler);
                if (tempDist == 0.) return 0.;
                dist = std::min(dist, tempDist);
            }
            return dist;
        },
        [&polygon, &ruler](const mapbox::geometry::line_string<double>& line) {
            return lineToLinesDistance(line, polygon, ruler);
        },
        [&polygon, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
            double dist = std::numeric_limits<double>::infinity();
            for (const auto& line : lines) {
                auto tempDist = lineToLinesDistance(line, polygon, ruler);
                if (tempDist == 0.) return 0.;
                dist = std::min(dist, tempDist);
            }
            return dist;
        },
        [&polygon, &ruler](const mapbox::geometry::polygon<double>& polygon1) {
            double dist = std::numeric_limits<double>::infinity();
            for (const auto& line : polygon1) {
                auto tempDist = lineToPolygonDistance(line, polygon, ruler);
                if (tempDist == 0.) return 0.;
                dist = std::min(dist, tempDist);
            }
            return dist;
        },
        [&polygon, &ruler](const mapbox::geometry::multi_polygon<double>& polygons) {
            double dist = std::numeric_limits<double>::infinity();
            for (const auto& poly : polygons) {
                for (const auto& line : poly) {
                    auto tempDist = lineToPolygonDistance(line, polygon, ruler);
                    if (tempDist == 0.) return 0.;
                    dist = std::min(dist, tempDist);
                }
            }
            return dist;
        },
        [](const auto&) -> double { return std::numeric_limits<double>::quiet_NaN(); });
}

double calculateDistance(const GeometryTileFeature& feature,
                         const CanonicalTileID& canonical,
                         const Feature::geometry_type& geoSet,
                         mapbox::cheap_ruler::CheapRuler::Unit unit) {
    return convertGeometry(feature, canonical)
        .match(
            [&geoSet, &unit](const mapbox::geometry::point<double>& point) -> double {
                return pointToGeometryDistance(point, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_point<double>& points) -> double {
                double ret = std::numeric_limits<double>::infinity();
                for (const auto& p : points) {
                    auto dist = pointToGeometryDistance(p, geoSet, unit);
                    if (dist == 0.) return dist;
                    ret = std::min(ret, dist);
                }
                return ret;
            },
            [&geoSet, &unit](const mapbox::geometry::line_string<double>& line) -> double {
                return lineToGeometryDistance(line, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_line_string<double>& lines) -> double {
                double ret = std::numeric_limits<double>::infinity();
                for (const auto& line : lines) {
                    auto dist = lineToGeometryDistance(line, geoSet, unit);
                    if (dist == 0.) return dist;
                    ret = std::min(ret, dist);
                }
                return ret;
            },
            [&geoSet, &unit](const mapbox::geometry::polygon<double>& polygon) -> double {
                return polygonToGeometryDistance(polygon, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_polygon<double>& polygons) -> double {
                double ret = std::numeric_limits<double>::infinity();
                for (const auto& polygon : polygons) {
                    auto dist = polygonToGeometryDistance(polygon, geoSet, unit);
                    if (dist == 0.) return dist;
                    ret = std::min(ret, dist);
                }
                return ret;
            },
            [](const auto&) -> double { return std::numeric_limits<double>::quiet_NaN(); });
}

struct Arguments {
    Arguments(GeoJSON& geojson_, mapbox::cheap_ruler::CheapRuler::Unit unit_)
        : geojson(std::move(geojson_)), unit(unit_) {}

    GeoJSON geojson;
    mapbox::cheap_ruler::CheapRuler::Unit unit;
};

optional<Arguments> parseValue(const style::conversion::Convertible& value, style::expression::ParsingContext& ctx) {
    if (isArray(value)) {
        // object value, quoted with ["Distance", GeoJSONObj, "unit"]
        auto length = arrayLength(value);
        if (length != 2 && length != 3) {
            ctx.error("'distance' expression requires exactly one argument, but found " +
                      util::toString(arrayLength(value) - 1) + " instead.");
            return nullopt;
        }

        mapbox::cheap_ruler::CheapRuler::Unit unit = mapbox::cheap_ruler::CheapRuler::Unit::Meters;
        if (arrayLength(value) == 3) {
            auto input = toString(arrayMember(value, 2)).value_or("Meters");
            if (input == "Meters" || input == "Metres") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Meters;
            }
            if (input == "Kilometers") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Kilometers;
            }
            if (input == "Miles") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Miles;
            }
            if (input == "Inches") {
                unit = mapbox::cheap_ruler::CheapRuler::Unit::Inches;
            }
        }
        const auto& argument1 = arrayMember(value, 1);
        if (isObject(argument1)) {
            style::conversion::Error error;
            auto geojson = toGeoJSON(argument1, error);
            if (geojson && error.message.empty()) {
                return Arguments(*geojson, unit);
            }
            ctx.error(error.message);
        }
    }
    ctx.error("'distance' expression needs to be an array with one/two arguments.");
    return nullopt;
}

optional<Feature::geometry_type> getGeometry(const Feature& feature, mbgl::style::expression::ParsingContext& ctx) {
    const auto type = apply_visitor(ToFeatureType(), feature.geometry);
    if (type == FeatureType::Point || type == FeatureType::LineString) {
        return feature.geometry;
    }
    ctx.error("'distance' expression requires valid geojson source that contains Point/LineString geometry type.");
    return nullopt;
}
} // namespace

namespace style {
namespace expression {

Distance::Distance(GeoJSON geojson, Feature::geometry_type geometries_, mapbox::cheap_ruler::CheapRuler::Unit unit_)
    : Expression(Kind::Distance, type::Number),
      geoJSONSource(std::move(geojson)),
      geometries(std::move(geometries_)),
      unit(unit_) {}

Distance::~Distance() = default;

using namespace mbgl::style::conversion;

EvaluationResult Distance::evaluate(const EvaluationContext& params) const {
    if (!params.feature || !params.canonical) {
        return EvaluationError{"distance expression requirs valid feature and canonical information."};
    }
    auto geometryType = params.feature->getType();
    if (geometryType == FeatureType::Unknown) return EvaluationError{"unknown geometry type."};
    auto distance = calculateDistance(*params.feature, *params.canonical, geometries, unit);
    if (std::isnan(distance)) {
        return EvaluationError{"unknown geometry type for evaluation."};
    }
    return distance;
}

ParseResult Distance::parse(const Convertible& value, ParsingContext& ctx) {
    auto parsedValue = parseValue(value, ctx);
    if (!parsedValue) {
        return ParseResult();
    }

    return parsedValue->geojson.match(
        [&parsedValue, &ctx](const mapbox::geometry::geometry<double>& geometrySet) {
            if (auto ret = getGeometry(mbgl::Feature(geometrySet), ctx)) {
                return ParseResult(
                    std::make_unique<Distance>(parsedValue->geojson, std::move(*ret), parsedValue->unit));
            }
            return ParseResult();
        },
        [&parsedValue, &ctx](const mapbox::feature::feature<double>& feature) {
            if (auto ret = getGeometry(mbgl::Feature(feature), ctx)) {
                return ParseResult(
                    std::make_unique<Distance>(parsedValue->geojson, std::move(*ret), parsedValue->unit));
            }
            return ParseResult();
        },
        [&parsedValue, &ctx](const mapbox::feature::feature_collection<double>& features) {
            for (const auto& feature : features) {
                if (auto ret = getGeometry(mbgl::Feature(feature), ctx)) {
                    return ParseResult(
                        std::make_unique<Distance>(parsedValue->geojson, std::move(*ret), parsedValue->unit));
                }
            }
            return ParseResult();
        },
        [&ctx](const auto&) {
            ctx.error("'distance' expression requires valid geojson that contains LineString/Point geometries.");
            return ParseResult();
        });

    return ParseResult();
}

Value convertValue(const mapbox::geojson::rapidjson_value& v) {
    if (v.IsDouble()) {
        return v.GetDouble();
    }
    if (v.IsString()) {
        return std::string(v.GetString());
    }
    if (v.IsArray()) {
        std::vector<Value> result;
        result.reserve(v.Size());
        for (const auto& m : v.GetArray()) {
            result.push_back(convertValue(m));
        }
        return result;
    }
    if (v.IsObject()) {
        std::unordered_map<std::string, Value> result;
        for (const auto& m : v.GetObject()) {
            result.emplace(m.name.GetString(), convertValue(m.value));
        }
        return result;
    }
    // Ignore other types as valid geojson only contains above types.
    return Null;
}

mbgl::Value Distance::serialize() const {
    std::unordered_map<std::string, Value> serialized;
    rapidjson::CrtAllocator allocator;
    const mapbox::geojson::rapidjson_value value = mapbox::geojson::convert(geoJSONSource, allocator);
    if (value.IsObject()) {
        for (const auto& m : value.GetObject()) {
            serialized.emplace(m.name.GetString(), convertValue(m.value));
        }
    } else {
        mbgl::Log::Error(mbgl::Event::General,
                         "Failed to serialize 'distance' expression, converted rapidJSON is not an object");
    }
    return std::vector<mbgl::Value>{{getOperator(), *fromExpressionValue<mbgl::Value>(serialized)}};
}

bool Distance::operator==(const Expression& e) const {
    if (e.getKind() == Kind::Distance) {
        auto rhs = static_cast<const Distance*>(&e);
        return geoJSONSource == rhs->geoJSONSource && geometries == rhs->geometries && unit == rhs->unit;
    }
    return false;
}

std::vector<optional<Value>> Distance::possibleOutputs() const {
    return {nullopt};
}

std::string Distance::getOperator() const {
    return "distance";
}

} // namespace expression
} // namespace style
} // namespace mbgl
