# 📊 SampleSizeCalculator

[![PyPI version](https://img.shields.io/pypi/v/sample-size-calculator.svg)](https://pypi.org/project/sample-size-calculator/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue)](https://www.python.org/)

`SampleSizeCalculator` — Python-класс для расчёта минимального размера выборки в A/B-тестах с помощью z-теста и поправки Бонферрони. Подходит для data-аналитиков, маркетологов и продакт-менеджеров.

---

## 🚀 Возможности

- Расчёт минимальной выборки с учётом:
  - Базовой конверсии
  - Минимального обнаруживаемого эффекта (MDE)
  - Баланса между тестом и контролем
  - Поправки на множественные сравнения
- Поддержка визуализации зависимости размера выборки от MDE
- Быстрая оценка доверительного интервала

---

## 📦 Установка

```bash
python3 -m pip install git+https://github.com/Aproject618/sample_size_calculator.git
```

Или, если устанавливаете из исходников:

```bash
git clone https://github.com/Aproject618/sample_size_calculator.git
cd sample-size-calculator
pip install .
```

## 🧪 Пример использования
```Python
from sample_size_calculator_packages import SampleSizeCalculator

# Инициализация калькулятора
calculator = SampleSizeCalculator(
    baseline_cr=0.05,    # Базовая конверсия (например, 5%)
    power=0.8,           # Мощность теста
    alpha=0.05,          # Уровень значимости
    ratio=0.5,           # Пропорция теста
    comparisons=1        # Количество гипотез (для Бонферрони)
)

# Расчёт размера выборки для 1% эффекта
sample_size, ci = calculator.calculate(mde=0.01)
print(f"Необходимый размер выборки: {int(sample_size)}")
print(f"Доверительный интервал: {ci}")
```

## 📈 Визуализация
```Python
calculator.plot_sample_sizes(mde_range=[5, 10, 15, 20, 25])
```
## 📋 Аргументы конструктора
| Аргумент      | Тип     | Описание                                  |
| ------------- | ------- | ----------------------------------------- |
| `baseline_cr` | `float` | Базовая конверсия (от 0 до 1)             |
| `power`       | `float` | Мощность теста (например, 0.8)            |
| `alpha`       | `float` | Уровень значимости (например, 0.05)       |
| `ratio`       | `float` | Пропорция трафика в тесте (например, 0.5) |
| `comparisons` | `int`   | Кол-во гипотез для поправки Бонферрони    |

## 📌 Применение

    Планирование A/B и сплит-тестов

    Оценка трафика, необходимого для выявления эффекта

    Учёт MDE и множественных гипотез в экспериментах

## 🔧 Зависимости

    scipy

    numpy

    plotly

## 🧾 Лицензия

Лицензировано под MIT License.