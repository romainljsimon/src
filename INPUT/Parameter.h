/*
MIT License
Copyright (c) 2019 - present H. Watanabe
The latest version is available at
https://github.com/kaityo256/params

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef PARAMETER_H_
#define PARAMETER_H_

#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace param {
    typedef std::map<std::string, std::string> pType;
    class Parameter {
    private:
        pType params;
        bool valid;

        bool contains(const std::string& key)
        {
            return params.find(key) != params.end();
        }

    public:
        explicit Parameter(const std::string& fileName)
                : valid(true)
        {
            loadFromFile(fileName);
        }

        explicit operator bool() const
        {
            return valid;
        };

        void loadFromFile(const std::string& filename)
        {
            std::ifstream is(filename);
            if (is.fail())
            {
                std::cerr << "Could not open file " << filename << std::endl;
                valid = false;
            }
            readFromStream(is);
        }
        void readFromStream(std::istream &is)
        {
            std::string line;
            while (getline(is, line))
            {
                if (line.length() > 0 && line[0] == '#')
                {
                    continue;
                }
                size_t index = line.find('=');
                if (std::string::npos != index)
                {
                    std::string key = line.substr(0, index);
                    std::string value = line.substr(index + 1, line.length());
                    params.insert(pType::value_type(key, value));
                }
            }
        }
        void check_key(const std::string& key) {
            if (!contains(key))
            {
                std::cerr << "No such key: " << key << "\n";
                std::abort();
            }
        }

        bool get_bool(const std::string& key, bool value)
        {
            if (!contains(key))
            {
                return value;
            }
            return ("yes" == params[key] || "Yes" == params[key]);
        }

        bool get_bool(const std::string& key)
        {
            check_key(key);
            return get_bool(key, false);
        }

        int get_int(const std::string& key, int value)
        {
            if (!contains(key))
            {
                return value;
            }
            return std::stoi(params[key]);
        }

        int get_int(const std::string& key)
        {
            check_key(key);
            return get_int(key, 0);
        }

        double get_double(const std::string& key, double value)
        {
            if (!contains(key)) {
                return value;
            }
            return std::stod(params[key]);
        }

        double get_double(const std::string& key)
        {
            check_key(key);
            return get_double(key, 0.0);
        }

        std::string get_string(const std::string& key, std::string value)
        {
            if (!contains(key))
            {
                return value;
            }
            return params[key];
        }

        std::string get_string(const std::string& key)
        {
            check_key(key);
            return get_string(key, "");
        }
    };



} // namespace param
#endif /* PARAMETER_H_ */
