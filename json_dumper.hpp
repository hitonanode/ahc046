#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <vector>

class JsonDumper {
    struct KeyValue {
        std::string key;
        std::string value;
    };

    std::vector<KeyValue> _items;

    bool dump_at_end = false;

public:
    JsonDumper(bool dump_at_end_ = false) : dump_at_end(dump_at_end_) {}

    ~JsonDumper() {
        if (dump_at_end) std::cout << dump() << std::endl;
    }

    void set_dump_at_end() { dump_at_end = true; }

    void operator()(const std::string &key, const std::string &value) {
        _items.push_back(KeyValue{key, "\"" + value + "\""});
    }

    template <class T> void operator()(const std::string &key, T value) {
        _items.push_back(KeyValue{key, std::to_string(value)});
    }

    std::string dump() const {
        std::string ret = "{\n";

        if (!_items.empty()) {
            for (const auto &[k, v] : _items) ret += "    \"" + k + "\": " + v + ",\n";

            ret.erase(ret.end() - 2);
        }

        ret += "}";
        return ret;
    }
} jdump;
