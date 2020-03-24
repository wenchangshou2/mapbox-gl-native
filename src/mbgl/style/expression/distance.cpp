#include <mbgl/style/expression/distance.hpp>

#include <mapbox/cheap_ruler.hpp>
#include <mapbox/geojson.hpp>
#include <mapbox/geometry.hpp>

#include <mbgl/style/conversion/json.hpp>
#include <mbgl/tile/geometry_tile_data.hpp>

#include <mbgl/util/logging.hpp>
#include <mbgl/util/string.hpp>

#include <rapidjson/document.h>

#include <tuple>

namespace mbgl {
namespace {

double pointDitanceToGeometry(const mapbox::geometry::point<double>& point,
                              const Feature::geometry_type& geoSet,
                              mapbox::cheap_ruler::CheapRuler::Unit unit) {
    mapbox::cheap_ruler::CheapRuler ruler(point.y, unit);
    return geoSet.match([&point, &ruler](const mapbox::geometry::point<double>& p) { return ruler.distance(point, p); },
                        [&point, &ruler](const mapbox::geometry::multi_point<double>& points) {
                            double ret = std::numeric_limits<double>::infinity();
                            for (size_t i = 0; i < points.size(); ++i) {
                                ret = std::min(ret, ruler.distance(point, points[i]));
                            }
                            return ret;
                        },
                        [&point, &ruler](const mapbox::geometry::line_string<double>& line) {
                            const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
                            return ruler.distance(point, nearestPoint);
                        },
                        [&point, &ruler](const mapbox::geometry::multi_line_string<double>& lines) {
                            double ret = std::numeric_limits<double>::infinity();
                            for (const auto& line : lines) {
                                const auto nearestPoint = std::get<0>(ruler.pointOnLine(line, point));
                                ret = std::min(ret, ruler.distance(point, nearestPoint));
                            }
                            return ret;
                        },
                        [](const auto&) -> double { return -1.0; });
}

double calculateDistance(const GeometryTileFeature& feature,
                         const CanonicalTileID& canonical,
                         const Feature::geometry_type& geoSet,
                         mapbox::cheap_ruler::CheapRuler::Unit unit) {
    return convertGeometry(feature, canonical)
        .match(
            [&geoSet, &unit](const mapbox::geometry::point<double>& point) -> double {
                return pointDitanceToGeometry(point, geoSet, unit);
            },
            [&geoSet, &unit](const mapbox::geometry::multi_point<double>& points) -> double {
                double ret = std::numeric_limits<double>::infinity();
                for (const auto& p : points) {
                    ret = std::min(ret, pointDitanceToGeometry(p, geoSet, unit));
                }
                return ret;
            },
            [](const auto&) -> double { return -1.0; });
}

struct Arguments {
    Arguments(GeoJSON& geojson_, mapbox::cheap_ruler::CheapRuler::Unit unit_)
        : geojson(std::move(geojson_)), unit(unit_) {}

    GeoJSON geojson;
    mapbox::cheap_ruler::CheapRuler::Unit unit;
};

optional<Arguments> parseValue(const style::conversion::Convertible& value, style::expression::ParsingContext& ctx) {
    if (isArray(value)) {
        // object value, quoted with ["Distance", GeoJSONObj, units]
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
        return -1.0;
    }
    auto geometryType = params.feature->getType();
    // Currently only support Point/Points distance to points or LineString
    if (geometryType == FeatureType::Point) {
        auto distance = calculateDistance(*params.feature, *params.canonical, geometries, unit);
        return distance;
    }
    mbgl::Log::Warning(mbgl::Event::General, "distance expression currently only support feature with Point geometry.");

    return -1.0;
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
